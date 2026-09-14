#!/usr/bin/env python3
"""Merge and deduplicate packing_graph.json datasets.

Graphs are deduplicated independently of run-local sample, atom, net, and node
IDs. Atom nodes are identified by their model and physical primitive number;
edges are canonicalized by endpoint and net grouping while preserving parallel
edge multiplicity.
"""

import argparse
import hashlib
import os
import sys
import tempfile
import time
from pathlib import Path
from typing import Any, BinaryIO, Callable, Iterator

try:
    import ijson
    import orjson
except ImportError as error:
    raise SystemExit(
        "deduplicate_packing_graphs.py requires ijson and orjson; install them "
        "with '.venv/bin/python -m pip install ijson orjson'"
    ) from error


READ_CHUNK_SIZE = 4 * 1024 * 1024


class PackingGraphError(RuntimeError):
    """Raised when an input packing graph dataset is malformed."""


def human_size(num_bytes: int) -> str:
    """Format a byte count for the progress display."""
    value = float(num_bytes)
    for unit in ("B", "KiB", "MiB", "GiB", "TiB"):
        if value < 1024.0 or unit == "TiB":
            return f"{value:.1f} {unit}"
        value /= 1024.0
    raise AssertionError("unreachable")


class ProgressBar:
    """Render throttled dataset progress to stderr."""

    def __init__(self, total_bytes: int, total_files: int, enabled: bool):
        self.total_bytes = total_bytes
        self.total_files = total_files
        self.enabled = enabled
        self.processed_bytes = 0
        self.last_render_time = 0.0

    def advance(
        self,
        num_bytes: int,
        completed_files: int,
        total_samples: int,
        unique_samples: int,
        duplicate_samples: int,
        force: bool = False,
    ) -> None:
        """Advance the byte count and redraw the bar when needed."""
        self.processed_bytes += num_bytes
        if not self.enabled:
            return

        current_time = time.monotonic()
        if not force and current_time - self.last_render_time < 0.1:
            return
        self.last_render_time = current_time

        fraction = (
            min(self.processed_bytes / self.total_bytes, 1.0)
            if self.total_bytes
            else 1.0
        )
        bar_width = 30
        filled_width = round(bar_width * fraction)
        bar = "#" * filled_width + "-" * (bar_width - filled_width)
        status = (
            f"\r[{bar}] {fraction:6.2%} "
            f"{human_size(self.processed_bytes)}/{human_size(self.total_bytes)} "
            f"files {completed_files}/{self.total_files} "
            f"samples {total_samples:,} unique {unique_samples:,} "
            f"duplicates {duplicate_samples:,}"
        )
        sys.stderr.write(status)
        sys.stderr.flush()

    def finish(self) -> None:
        """Terminate the progress line."""
        if self.enabled:
            sys.stderr.write("\n")
            sys.stderr.flush()


class ProgressReader:
    """Report byte progress while presenting a binary stream to ijson."""

    def __init__(self, input_file: BinaryIO, callback: Callable[[int], None]):
        self.input_file = input_file
        self.callback = callback

    def read(self, size: int = -1) -> bytes:
        """Read bytes and report the amount consumed."""
        chunk = self.input_file.read(size)
        self.callback(len(chunk))
        return chunk


def iter_packing_graphs(
    path: Path, progress_callback: Callable[[int], None] | None = None
) -> Iterator[dict[str, Any]]:
    """Yield graph records from a packing_graph.json file without loading it all."""
    callback = progress_callback or (lambda _: None)
    with path.open("rb") as input_file:
        reader = ProgressReader(input_file, callback)
        try:
            values = ijson.items(
                reader, "item", use_float=True, buf_size=READ_CHUNK_SIZE
            )
            for index, value in enumerate(values):
                if not isinstance(value, dict):
                    raise PackingGraphError(
                        f"{path}: sample {index} is not a JSON object"
                    )
                yield value
        except ijson.JSONError as error:
            raise PackingGraphError(f"{path}: invalid JSON: {error}") from error


def stable_json(value: Any) -> bytes:
    """Return a deterministic compact JSON representation."""
    return orjson.dumps(value, option=orjson.OPT_SORT_KEYS)


def canonical_graph_digest(record: dict[str, Any]) -> bytes:
    """Return an ID-independent digest of a packing graph's input features."""
    try:
        nodes = record["nodes"]
        edges = record["edges"]
    except KeyError as error:
        raise PackingGraphError(f"record is missing field {error}") from error

    if not isinstance(nodes, list) or not isinstance(edges, list):
        raise PackingGraphError("record fields 'nodes' and 'edges' must be arrays")

    external_node_ids: dict[str, int] = {}
    atom_nodes: list[tuple[bytes, int, dict[str, Any]]] = []
    seen_node_ids: set[int] = set()

    for node in nodes:
        if not isinstance(node, dict) or "id" not in node or "type" not in node:
            raise PackingGraphError("every node must contain 'id' and 'type'")

        node_id = node["id"]
        node_type = node["type"]
        if not isinstance(node_id, int) or node_id in seen_node_ids:
            raise PackingGraphError("node IDs must be unique integers")
        seen_node_ids.add(node_id)

        attributes = {
            key: value
            for key, value in node.items()
            if key not in ("id", "atom_id")
        }
        attribute_key = stable_json(attributes)

        if node_type == "atom":
            atom_nodes.append((attribute_key, node_id, attributes))
        elif node_type in ("external_input", "external_output"):
            if node_type in external_node_ids:
                raise PackingGraphError(f"multiple {node_type!r} nodes found")
            external_node_ids[node_type] = node_id
        else:
            raise PackingGraphError(f"unsupported node type {node_type!r}")

    if set(external_node_ids) != {"external_input", "external_output"}:
        raise PackingGraphError(
            "each record must have one external_input and one external_output node"
        )

    atom_nodes.sort()
    for previous, current in zip(atom_nodes, atom_nodes[1:]):
        if previous[0] == current[0]:
            raise PackingGraphError(
                "cannot canonicalize atoms with identical model and primitive features"
            )

    canonical_node_ids = {
        external_node_ids["external_input"]: 0,
        external_node_ids["external_output"]: 1,
    }
    for canonical_id, (_, original_id, _) in enumerate(atom_nodes, start=2):
        canonical_node_ids[original_id] = canonical_id

    net_edges: dict[Any, list[tuple[bytes, dict[str, Any]]]] = {}
    for edge in edges:
        if not isinstance(edge, dict):
            raise PackingGraphError("every edge must be a JSON object")
        try:
            source = canonical_node_ids[edge["src"]]
            sink = canonical_node_ids[edge["dst"]]
            net_key = edge["atom_net_id"]
        except KeyError as error:
            raise PackingGraphError(
                f"edge references a missing field or node {error}"
            ) from error

        attributes = {
            key: value
            for key, value in edge.items()
            if key not in ("src", "dst", "atom_net_id", "external_atom_id")
        }
        attributes["src"] = source
        attributes["dst"] = sink
        net_edges.setdefault(net_key, []).append(
            (stable_json(attributes), attributes)
        )

    keyed_nets: list[tuple[bytes, list[dict[str, Any]]]] = []
    for group in net_edges.values():
        canonical_edges = [
            attributes for _, attributes in sorted(group, key=lambda item: item[0])
        ]
        keyed_nets.append((stable_json(canonical_edges), canonical_edges))
    canonical_nets = [
        net for _, net in sorted(keyed_nets, key=lambda item: item[0])
    ]
    graph_metadata = {
        key: value
        for key, value in record.items()
        if key not in ("sample_id", "result", "nodes", "edges")
    }
    canonical_form = {
        "metadata": graph_metadata,
        "atom_nodes": [attributes for _, _, attributes in atom_nodes],
        "nets": canonical_nets,
    }
    return hashlib.sha256(stable_json(canonical_form)).digest()


def parse_args() -> argparse.Namespace:
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "input_folder",
        type=Path,
        help="folder searched recursively for packing graph JSON files",
    )
    parser.add_argument(
        "output_file",
        nargs="?",
        type=Path,
        default=Path("packing_graph_deduplicated.json"),
        help="merged output file (default: packing_graph_deduplicated.json)",
    )
    parser.add_argument(
        "--pattern",
        default="packing_graph.json",
        help="recursive input filename pattern (default: packing_graph.json)",
    )
    parser.add_argument(
        "--keep-label-conflicts",
        action="store_true",
        help="keep one sample for each distinct label when equivalent graphs conflict",
    )
    parser.add_argument(
        "--no-progress",
        action="store_true",
        help="disable the interactive progress bar",
    )
    return parser.parse_args()


def deduplicate(args: argparse.Namespace) -> None:
    """Merge input datasets while removing equivalent graph samples."""
    input_folder = args.input_folder.resolve()
    output_file = args.output_file.resolve()
    if not input_folder.is_dir():
        raise PackingGraphError(f"input folder does not exist: {input_folder}")

    input_files = sorted(input_folder.rglob(args.pattern))
    if not input_files:
        raise PackingGraphError(
            f"no files matching {args.pattern!r} found under {input_folder}"
        )

    output_file.parent.mkdir(parents=True, exist_ok=True)
    seen_labels: dict[bytes, set[bytes]] = {}
    total_samples = 0
    unique_samples = 0
    duplicate_samples = 0
    label_conflicts = 0
    total_bytes = sum(path.stat().st_size for path in input_files)
    progress = ProgressBar(
        total_bytes,
        len(input_files),
        enabled=sys.stderr.isatty() and not args.no_progress,
    )

    temporary_path: Path | None = None
    try:
        with tempfile.NamedTemporaryFile(
            mode="wb",
            dir=output_file.parent,
            prefix=f".{output_file.name}.",
            suffix=".tmp",
            delete=False,
        ) as output:
            temporary_path = Path(output.name)
            output.write(b"[\n")

            for file_index, input_file in enumerate(input_files):
                file_samples = 0
                file_bytes_read = 0

                def update_file_progress(num_bytes: int) -> None:
                    nonlocal file_bytes_read
                    file_bytes_read += num_bytes
                    progress.advance(
                        num_bytes,
                        completed_files=file_index,
                        total_samples=total_samples,
                        unique_samples=unique_samples,
                        duplicate_samples=duplicate_samples,
                    )

                for record in iter_packing_graphs(input_file, update_file_progress):
                    file_samples += 1
                    total_samples += 1
                    if "result" not in record:
                        raise PackingGraphError(
                            f"{input_file}: sample {file_samples - 1} has no result"
                        )
                    graph_digest = canonical_graph_digest(record)
                    label = stable_json(record["result"])
                    labels = seen_labels.get(graph_digest)

                    if labels is not None and label in labels:
                        duplicate_samples += 1
                        continue

                    if labels is not None:
                        label_conflicts += 1
                        if not args.keep_label_conflicts:
                            raise PackingGraphError(
                                "equivalent graphs have conflicting results; "
                                f"conflict encountered in {input_file}, sample "
                                f"{record.get('sample_id', file_samples - 1)}"
                            )
                        labels.add(label)
                    else:
                        seen_labels[graph_digest] = {label}

                    output_record = dict(record)
                    output_record["sample_id"] = unique_samples
                    if unique_samples != 0:
                        output.write(b",\n")
                    output.write(orjson.dumps(output_record))
                    unique_samples += 1

                    if total_samples % 1000 == 0:
                        progress.advance(
                            0,
                            completed_files=file_index,
                            total_samples=total_samples,
                            unique_samples=unique_samples,
                            duplicate_samples=duplicate_samples,
                        )

                file_size = input_file.stat().st_size
                progress.advance(
                    file_size - file_bytes_read,
                    completed_files=file_index + 1,
                    total_samples=total_samples,
                    unique_samples=unique_samples,
                    duplicate_samples=duplicate_samples,
                    force=True,
                )
                if not progress.enabled:
                    print(
                        f"Read {file_samples:,} samples from {input_file}",
                        flush=True,
                    )

            output.write(b"\n]\n")

        os.replace(temporary_path, output_file)
        temporary_path = None
    finally:
        progress.finish()
        if temporary_path is not None:
            temporary_path.unlink(missing_ok=True)

    print(f"Input files:       {len(input_files):,}")
    print(f"Input samples:     {total_samples:,}")
    print(f"Unique samples:    {unique_samples:,}")
    print(f"Duplicates removed:{duplicate_samples:>12,}")
    print(f"Label conflicts:   {label_conflicts:>12,}")
    print(f"Output: {output_file}")


def main() -> int:
    """Run the command-line tool."""
    try:
        deduplicate(parse_args())
    except (OSError, PackingGraphError, TypeError, ValueError) as error:
        print(f"error: {error}", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
