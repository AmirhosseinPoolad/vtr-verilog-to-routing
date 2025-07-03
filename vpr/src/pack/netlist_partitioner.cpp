#include "netlist_partitioner.h"
#include "vpr_error.h"
#include "pack_types.h"
#include "prepack.h"
#include "flat_placement_types.h"
#include "vpr_context.h"
#include "vtr_assert.h"

// If there are less than 'partition_threshold' number of atoms of a model type, they are not partitioned.
// This is done to avoid QoR losses for hard blocks such as RAM and DSP blocks because of partitioning.
static constexpr int partition_threshold = 2500;

NetlistPartitioner::NetlistPartitioner(const FlatPlacementInfo& flat_placement_info, const Prepacker& prepacker, const AtomContext& atom_context, const DeviceContext& device_context)
    : prepacker_(prepacker)
    , flat_placement_info_(flat_placement_info)
    , atom_context_(atom_context)
    , device_context_(device_context) {
        for (AtomBlockId atom_blk : atom_context_.netlist().blocks()) {
            char* model_name = device_context_.arch->models.get_model(atom_context_.netlist().block_model(atom_blk)).name;
            model_count_[model_name]++;
        }
    }

NetlistPartition NetlistPartitioner::get_netlist_partition(e_partition_type partition_type, int num_partitions) {
    NetlistPartition partition_map(num_partitions);
    switch (partition_type)
    {
    case e_partition_type::SPATIAL:
        partition_map = get_spatial_partitioning(num_partitions);
        break;

    case e_partition_type::NONE:
        partition_map = get_unity_partitioning();
        break;
    
    default:
        VPR_FATAL_ERROR(VPR_ERROR_PACK, "Unknown netlist partition type selected: %d\n", (int)partition_type);
    }
    return partition_map;
}

bool NetlistPartitioner::should_partition_mol(PackMoleculeId mol_id) {
    for (AtomBlockId blk_id : prepacker_.get_molecule(mol_id).atom_block_ids) {
        if (!blk_id.is_valid()) {
            continue;
        }
        
        char* model_name = device_context_.arch->models.get_model(atom_context_.netlist().block_model(blk_id)).name;
        
        if(prepacker_.get_expected_lowest_cost_pb_gnode(blk_id)->pb_type->class_type == MEMORY_CLASS) {
            return false;
        }

        if (model_count_[model_name] >= partition_threshold) {
            return true;
        }
    }
    return false;
}

static std::pair<int, int> get_closest_factors(int num) {
    int sqrt = std::sqrt(num);
    while (num % sqrt != 0) {
        sqrt--;
    }
    return std::make_pair(sqrt, num / sqrt);
}

static std::pair<t_flat_pl_loc, t_flat_pl_loc> get_flat_placement_bounding_box(const FlatPlacementInfo& flat_placement_info, const AtomContext& atom_context) {
    constexpr float MAX_FLOAT = std::numeric_limits<float>::max();
    constexpr float MIN_FLOAT = std::numeric_limits<float>::lowest();
    t_flat_pl_loc min_coords({MAX_FLOAT, MAX_FLOAT, MAX_FLOAT});
    t_flat_pl_loc max_coords({MIN_FLOAT, MIN_FLOAT, MIN_FLOAT});

    for (AtomBlockId blk_id : atom_context.netlist().blocks()) {
        auto cur_blk_loc = flat_placement_info.get_pos(blk_id);
        max_coords.x = (cur_blk_loc.x > max_coords.x) ? cur_blk_loc.x : max_coords.x;
        max_coords.y = (cur_blk_loc.y > max_coords.y) ? cur_blk_loc.y : max_coords.y;
        max_coords.layer = (cur_blk_loc.layer > max_coords.layer) ? cur_blk_loc.layer : max_coords.layer;
        min_coords.x = (cur_blk_loc.x < min_coords.x) ? cur_blk_loc.x : min_coords.x;
        min_coords.y = (cur_blk_loc.y < min_coords.y) ? cur_blk_loc.y : min_coords.y;
        min_coords.layer = (cur_blk_loc.layer < min_coords.layer) ? cur_blk_loc.layer : min_coords.layer;
    }
    return std::make_pair(min_coords, max_coords);
}

NetlistPartition NetlistPartitioner::get_spatial_partitioning(int num_partitions) {
    auto [min_coords, max_coords] = get_flat_placement_bounding_box(flat_placement_info_, atom_context_);
    NetlistPartition partition_map(num_partitions);

    if (num_partitions == 1) {
        // No need to partition, return a unity partition
        return get_unity_partitioning();
    }

    auto axis_partitions = get_closest_factors(num_partitions);

    int num_partitions_x = axis_partitions.first;
    int num_partitions_y = axis_partitions.second;
    int num_partitions_layer = 1;

    double partition_size_x = (max_coords.x - min_coords.x) / num_partitions_x;
    double partition_size_y = (max_coords.y - min_coords.y) / num_partitions_y;
    int partition_size_layer = 1;
    VTR_LOG("Netlist Partition size: (%lf, %lf, %d)\n", partition_size_x, partition_size_y, partition_size_layer);

    // Get the partition for each molecule
    for (auto mol : prepacker_.molecules()) {
        if (!should_partition_mol(mol)) {
            partition_map.set_molecule_partition(mol, 0);
            continue;
        }
        auto cur_blk_loc = flat_placement_info_.get_pos(prepacker_.get_molecule(mol).atom_block_ids[0]);

        auto x = cur_blk_loc.x;
        auto y = cur_blk_loc.y;

        int partition_x = (int)((x - min_coords.x) / partition_size_x);
        int partition_y = (int)((y - min_coords.y) / partition_size_y);
        partition_x = partition_x == num_partitions_x ? partition_x - 1 : partition_x;
        partition_y = partition_y == num_partitions_y ? partition_y - 1 : partition_y;

        int partition_layer = 0;

        int partition_id = partition_x + (partition_y * num_partitions_x);

        // Check if the partition is valid. This should never happen
        if (partition_x < 0 || partition_x >= num_partitions_x || partition_y < 0 || partition_y >= num_partitions_y || partition_layer < 0 || partition_layer >= num_partitions_layer || partition_id < 0 || partition_id >= num_partitions) {
            VPR_FATAL_ERROR(VPR_ERROR_PACK, "Partition: (%d, %d, %d) = (%d) is out of bounds: max partitions(%d, %d, %d) = (%d)\n", partition_x, partition_y, partition_layer, partition_id, num_partitions_x, num_partitions_y, num_partitions_layer, num_partitions);
        }
        partition_map.set_molecule_partition(mol, partition_id);
    }
    return partition_map;
}

NetlistPartition NetlistPartitioner::get_unity_partitioning() {
    NetlistPartition partition_map(1);
    for (auto mol : prepacker_.molecules()) {
        partition_map.set_molecule_partition(mol, 0);
    }
    return partition_map;
}

const std::vector<PackMoleculeId>& NetlistPartition::molecules(int partition) const {
    VTR_ASSERT_SAFE((size_t)partition < molecules_.size());
    return molecules_[partition];
}

int NetlistPartition::get_partition(PackMoleculeId molecule) const {
    return partition_map_.at(molecule);
}

void NetlistPartition::set_molecule_partition(PackMoleculeId mol, int partition) {
    VTR_ASSERT_SAFE((size_t)partition < molecules_.size());
    if(partition_map_.contains(mol)) {
        int previous_partition = partition_map_[mol];
        std::vector<PackMoleculeId>& prev_partition_vec = molecules_[previous_partition];
        prev_partition_vec.erase(std::remove(prev_partition_vec.begin(), prev_partition_vec.end(), mol), prev_partition_vec.end());
    }
    partition_map_[mol] = partition;
    molecules_[partition].push_back(mol);
}

NetlistPartition::NetlistPartition(int num_partitions) {
    molecules_.resize(num_partitions);
}
