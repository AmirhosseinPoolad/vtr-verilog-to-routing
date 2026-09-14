#!/usr/bin/env bash

set -euo pipefail

readonly DEFAULT_SOURCE="/home/amirpoolad/Dev/vtr-verilog-to-routing/vtr_flow/tasks/regression_tests/vtr_reg_nightly_test2/vtr_xilinx_qor/run001"
readonly DEFAULT_DESTINATION="/home/amirpoolad/Dev/vtr-verilog-to-routing/graph_data"

source_directory="${1:-$DEFAULT_SOURCE}"
destination_directory="${2:-$DEFAULT_DESTINATION}"

if [[ ! -d "$source_directory" ]]; then
    echo "error: source directory does not exist: $source_directory" >&2
    exit 1
fi

mkdir -p "$destination_directory"

copied_count=0
while IFS= read -r -d '' source_file; do
    relative_path="${source_file#"$source_directory"/}"
    destination_file="$destination_directory/$relative_path"

    mkdir -p "$(dirname "$destination_file")"
    cp -- "$source_file" "$destination_file"
    ((copied_count += 1))
done < <(find "$source_directory" -type f -name packing_graph.json -print0)

echo "Copied $copied_count packing_graph.json file(s) to $destination_directory"
