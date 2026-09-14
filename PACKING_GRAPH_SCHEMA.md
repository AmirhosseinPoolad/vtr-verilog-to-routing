# Packing Graph JSON Schema

`packing_graph.json` contains labeled atom-level cluster graphs collected during
VPR packing. The dataset is intended for graph machine-learning tasks that
predict whether an intra-cluster routing attempt succeeds.

## File structure

The top-level value is a JSON array. Each element is one independent graph
sample:

```json
[
  {
    "schema": "vtr.intra_lb_route.v1",
    "sample_id": 0,
    "cluster_type": "CLB",
    "result": {
      "is_impossible": false,
      "route_succeeded": true
    },
    "nodes": [
      {"id": 0, "type": "external_input"},
      {"id": 1, "type": "external_output"},
      {"id": 2, "type": "atom", "atom_id": 31, "model": ".names", "primitive_num": 4}
    ],
    "edges": [
      {"src": 0, "dst": 2, "atom_net_id": 7},
      {"src": 2, "dst": 1, "atom_net_id": 8}
    ]
  }
]
```

## Sample fields

| Field | Type | Description |
| --- | --- | --- |
| `schema` | string | Schema identifier. The current value is `vtr.intra_lb_route.v1`. |
| `sample_id` | integer | Zero-based record number within this JSON file. It is not globally stable. |
| `cluster_type` | string | Name of the logical block type whose internal routing is being attempted. |
| `result` | object | Labels produced by the intra-cluster router. |
| `nodes` | array | Nodes in the directed graph. |
| `edges` | array | Directed atom-to-atom or boundary-to-atom connections. |

### Result labels

| Field | Type | Description |
| --- | --- | --- |
| `is_impossible` | boolean | `true` when at least one required connection has no legal path through the cluster routing-resource graph. |
| `route_succeeded` | boolean | `true` when the router succeeds on its final iteration. |

Samples rejected because of a mode-selection issue are not written. For the
remaining samples, `is_impossible` distinguishes structurally impossible
instances from possible instances. `route_succeeded` supplies the final routing
outcome; an impossible sample is expected to have `route_succeeded: false`.

## Nodes

Every sample contains exactly two boundary nodes followed by zero or more atom
nodes. Node IDs are dense and local to the sample.

### Boundary nodes

```json
{"id": 0, "type": "external_input"}
{"id": 1, "type": "external_output"}
```

- Node `0` represents all drivers outside the cluster. An edge from node `0`
  to an atom means that the atom consumes a net driven externally.
- Node `1` represents all sinks outside the cluster. An edge from an atom to
  node `1` means that the atom drives at least one sink outside the cluster.

The boundary nodes aggregate external atoms; they do not identify individual
external drivers or sinks.

### Atom nodes

```json
{
  "id": 2,
  "type": "atom",
  "atom_id": 31,
  "model": ".names",
  "primitive_num": 4
}
```

| Field | Type | Description |
| --- | --- | --- |
| `id` | integer | Dense node ID used by edge endpoints. Atom IDs begin at `2`. |
| `type` | string | Always `atom` for an atom node. |
| `atom_id` | integer | Atom block ID from the current VPR atom netlist. This is run-local metadata. |
| `model` | string | Logical model name, such as `.names`, `.latch`, or an architecture model. |
| `primitive_num` | integer | Cluster-relative physical primitive number assigned to the atom. A value of `-1` means no primitive number was recovered. |

## Edges

Each edge is directed from a net driver to one sink:

```json
{"src": 2, "dst": 3, "atom_net_id": 42}
```

| Field | Type | Description |
| --- | --- | --- |
| `src` | integer | Source node ID. Must refer to a node in the same sample. |
| `dst` | integer | Destination node ID. Must refer to a node in the same sample. |
| `atom_net_id` | integer | Atom netlist net ID. This is run-local metadata. Equal values identify edges expanded from the same hyperedge. |

An atom-net hyperedge is expanded into one directed edge per sink. Therefore,
multiple edges may share the same source and `atom_net_id`. Parallel edges are
valid and their multiplicity is meaningful; consumers should not collapse them
unless that behavior is explicitly desired.

The supported endpoint forms are:

| Connection | `src` | `dst` |
| --- | --- | --- |
| Internal connection | atom | atom |
| External input | `external_input` | atom |
| External output | atom | `external_output` |

## Identifier scope

`sample_id`, `atom_id`, `atom_net_id`, and node `id` are identifiers rather
than learning features:

- `sample_id` is meaningful only within one file.
- `atom_id` and `atom_net_id` come from one VPR invocation and may change
  between equivalent runs.
- Node `id` is meaningful only inside its sample and is referenced by `src` and
  `dst`.
- `primitive_num` is a graph feature: it identifies placement within the
  cluster's physical primitive structure.

## Deduplicated datasets

`deduplicate_packing_graphs.py` accepts a directory tree containing files with
this schema and writes one top-level JSON array with the same record structure.
It replaces `sample_id` with a new dense sequence in the merged output.

Deduplication ignores run-local sample, atom, net, and node identifiers. It
preserves cluster metadata, atom features, directed topology, net grouping,
parallel-edge multiplicity, and result labels. By default, equivalent graphs
with conflicting result labels are reported as an error; the
`--keep-label-conflicts` option retains one copy of each distinct label.
