#include "netlist_partitioner.h"
#include "vpr_error.h"
#include "log.h"
#include "pack_types.h"
#include "prepack.h"
#include "flat_placement_types.h"
#include "vpr_context.h"
#include "vtr_assert.h"
#include <mtkahypar.h>

#include <iostream>
#include <fstream>

#include <sys/types.h>
#include <sys/wait.h>
#include <unistd.h>
#include <iostream>

static constexpr int partition_threshold = 5000;

static std::pair<int, int> get_closest_factors(int num) {
    int sqrt = std::sqrt(num);
    while (num % sqrt != 0) {
        sqrt--;
    }
    return std::make_pair(sqrt, num / sqrt);
}

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
    NetlistPartition partition_map(t_partition_dimension(0,0,0));
    switch (partition_type)
    {
    case e_partition_type::MIN_CUT:
        partition_map = get_graph_partitioning(num_partitions, false);
        break;

    case e_partition_type::SPATIAL_MIN_CUT:
        partition_map = get_graph_partitioning(num_partitions, true);
        break;

    case e_partition_type::SPATIAL:
        partition_map = get_spatial_partitioning(num_partitions);
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

        if(prepacker_.get_expected_lowest_cost_pb_gnode(blk_id)->pb_type->class_type == MEMORY_CLASS) {
            return false;
        }

        char* model_name = device_context_.arch->models.get_model(atom_context_.netlist().block_model(blk_id)).name;
        // TODO: magic number
        if (model_count_[model_name] >= partition_threshold) {
            return true;
        }
    }
    return false;
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

    // VTR_LOG("Flat placement min coords: (%g, %g, %g)\n", min_coords.x, min_coords.y, min_coords.layer);
    // VTR_LOG("Flat placement max coords: (%g, %g, %g)\n", max_coords.x, max_coords.y, max_coords.layer);
    // VTR_LOG("Flat placement size: (%g, %g, %g)\n", max_coords.x - min_coords.x, max_coords.y - min_coords.y, max_coords.layer - min_coords.layer);

    return std::make_pair(min_coords, max_coords);
}

NetlistPartition NetlistPartitioner::get_spatial_partitioning(int num_partitions) {
    auto [min_coords, max_coords] = get_flat_placement_bounding_box(flat_placement_info_, atom_context_);

    auto axis_partitions = get_closest_factors(num_partitions);

    int num_partitions_x = axis_partitions.first;
    int num_partitions_y = axis_partitions.second;
    int num_partitions_z = 1;

    t_partition_dimension partition_dimensions = {.x = num_partitions_x, .y = num_partitions_y, .z = num_partitions_z};

    NetlistPartition partition_map(partition_dimensions);

    if (num_partitions == 1) {
        // No need to partition, just return 0
        for (auto mol : prepacker_.molecules()) {
            partition_map.set_molecule_partition(mol, 0);
        }
        return partition_map;
    }

    

    double partition_size_x = (max_coords.x - min_coords.x) / partition_dimensions.x;
    double partition_size_y = (max_coords.y - min_coords.y) / partition_dimensions.y;
    int partition_size_layer = 1;
    VTR_LOG("Partition size: (%lf, %lf, %d)\n", partition_size_x, partition_size_y, partition_size_layer);

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
        partition_x = partition_x == partition_dimensions.x ? partition_x - 1 : partition_x;
        partition_y = partition_y == partition_dimensions.y ? partition_y - 1 : partition_y;

        int partition_layer = 0;

        int partition_id = partition_x + (partition_y * partition_dimensions.x);

        // Check if the partition is valid. This should never happen
        if (partition_x < 0 || partition_x >= partition_dimensions.x || partition_y < 0 || partition_y >= partition_dimensions.y || partition_layer < 0 || partition_layer >= partition_dimensions.z || partition_id < 0 || partition_id >= num_partitions) {
            VPR_FATAL_ERROR(VPR_ERROR_PACK, "Partition: (%d, %d, %d) = (%d) is out of bounds: max partitions(%d, %d, %d) = (%d)\n", partition_x, partition_y, partition_layer, partition_id, partition_dimensions.x, partition_dimensions.y, partition_dimensions.z, num_partitions);
        }
        partition_map.set_molecule_partition(mol, partition_id);
    }
    return partition_map;
}

static float compute_mean_distance(const std::unordered_map<PackMoleculeId, std::pair<float, float>>& block_positions) {
    float total_distance = 0.0f;
    int count = 0;

    for (auto it1 = block_positions.begin(); it1 != block_positions.end(); ++it1) {
        auto it2 = it1;
        ++it2;
        for (; it2 != block_positions.end(); ++it2) {
            float dx = it1->second.first - it2->second.first;
            float dy = it1->second.second - it2->second.second;
            float dist = std::sqrt(dx * dx + dy * dy);

            total_distance += dist;
            ++count;
        }
    }

    return (count > 0) ? (total_distance / count) : 0.0f;
}

static void write_hmetis_file(bool use_placement_info, const AtomContext& atom_context, const Prepacker& prepacker, const FlatPlacementInfo& flat_placement_info) {
    std::ofstream graph_traverse_file("graph_traverse.txt", std::ios::out);
    if (!use_placement_info) {
        graph_traverse_file << atom_context.netlist().nets().size() << " " << prepacker.molecules().size() << std::endl;
    }
    std::vector<AtomBlockId> block_ids;
    std::vector<PackMoleculeId> molecule_ids;
    std::unordered_map<PackMoleculeId, std::pair<float, float>> block_positions;
    int inter_molecule_hyperedges = 0;
    // Create graph of molecules
    bool put_space = false;
    if (graph_traverse_file.is_open()) {
        for (auto net_id : atom_context.netlist().nets()) {
            molecule_ids.clear();
            put_space = false;
            for (auto pin_id : atom_context.netlist().net_pins(net_id)) {
                auto block_id = atom_context.netlist().pin_block(pin_id);
                auto target_molecule_id = prepacker.get_atom_molecule(block_id);
                auto it = std::find(molecule_ids.begin(), molecule_ids.end(), target_molecule_id);
                if (it == molecule_ids.end()) {
                    if (use_placement_info) {
                        if (block_positions.find(target_molecule_id) == block_positions.end()) {
                            auto molecule_root_block_id = prepacker.get_molecule_root_atom(target_molecule_id);
                            // Get the position of the root block in this molecule
                            auto [block_x, block_y, block_layer] = flat_placement_info.get_pos(molecule_root_block_id);
                            block_positions[target_molecule_id] = {block_x, block_y};
                        }
                    }
                    molecule_ids.push_back(target_molecule_id);
                }
            }

            if (molecule_ids.size() > 1) {
                if (use_placement_info) {
                    inter_molecule_hyperedges++;
                    float hyperedge_weight = compute_mean_distance(block_positions);
                    constexpr float distance_weight_coeff = 10.0f;
                    graph_traverse_file << (int)(distance_weight_coeff/hyperedge_weight);
                    for (auto target_molecule_id : molecule_ids) {
                        graph_traverse_file << " ";                    
                        graph_traverse_file << ((int)target_molecule_id + 1);
                        }   
                        graph_traverse_file << "\n"; // Only write inter-molecule hyperedges
                    }                         
                
                else {
                    for (auto target_molecule_id : molecule_ids) {
                        if (put_space) {
                            graph_traverse_file << " ";  
                        }                    
                            graph_traverse_file << ((int)target_molecule_id + 1);
                            put_space = true;
                        }       
                    } 
            }                  
            if (!use_placement_info) {
                graph_traverse_file << "\n";
            }
            block_positions.clear();       
        }
        graph_traverse_file.close();
        if (use_placement_info) {
            std::ifstream in_file("graph_traverse.txt");
            std::ostringstream buffer;
            buffer << in_file.rdbuf();
            in_file.close();     
            std::string content = buffer.str();
            if (!content.empty() && content.back() == '\n') {
                content.pop_back(); // remove only the last newline
            }
            
            // FIXME: why is this declared here again
            std::ofstream graph_traverse_file2("graph_traverse.txt", std::ios::out);
            graph_traverse_file2 << inter_molecule_hyperedges << " " << prepacker.molecules().size() << " " << "1" << std::endl;
            graph_traverse_file2 << content;
            graph_traverse_file2.close();
        }
    } else {
        VPR_FATAL_ERROR(VPR_ERROR_PACK, "Error opening hmetis file\n");
    }
}

// mt_kahypar is for some reason broken and will corrupt the heap
// call this function in a seperate process
// FIXME: Use the compiled binary instead of the library
static void call_mt_kahypar(int num_partitions) {
    mt_kahypar_error_t error{};
    // Initialize MT-KaHyPar
    mt_kahypar_initialize(1, true);
    mt_kahypar_context_t* context = mt_kahypar_context_from_preset(DEFAULT);
    mt_kahypar_set_partitioning_parameters(context, num_partitions, 0.10, CUT);
    mt_kahypar_set_seed(42);
    mt_kahypar_status_t status = mt_kahypar_set_context_parameter(context, VERBOSE, "1", &error);
    VTR_ASSERT(status == SUCCESS);

    mt_kahypar_hypergraph_t hypergraph = mt_kahypar_read_hypergraph_from_file("graph_traverse.txt", context, HMETIS, &error);
    if (error.status != SUCCESS || hypergraph.hypergraph == nullptr) {
        std::cerr << error.msg << std::endl;
        std::exit(1);
    }

    mt_kahypar_partitioned_hypergraph_t partitioned_hg = mt_kahypar_partition(hypergraph, context, &error);
    if (partitioned_hg.partitioned_hg == nullptr) {
        std::cerr << error.msg << std::endl;
        std::exit(1);
    }

    auto partition = std::make_unique<mt_kahypar_partition_id_t[]>(mt_kahypar_num_hypernodes(hypergraph));
    mt_kahypar_get_partition(partitioned_hg, partition.get());

    auto block_weights = std::make_unique<mt_kahypar_hypernode_weight_t[]>(2);
    mt_kahypar_get_block_weights(partitioned_hg, block_weights.get());

    const double imbalance = mt_kahypar_imbalance(partitioned_hg, context);
    const int km1 = mt_kahypar_km1(partitioned_hg);
    mt_kahypar_write_partition_to_file(partitioned_hg, "partitioned_graph.txt", &error);

    VTR_ASSERT(error.status == SUCCESS);

    std::cout << "Partitioning Results:" << std::endl;
    std::cout << "Imbalance         = " << imbalance << std::endl;
    std::cout << "Km1               = " << km1 << std::endl;
    std::cout << "Weight of Block 0 = " << block_weights[0] << std::endl;
    std::cout << "Weight of Block 1 = " << block_weights[1] << std::endl;

    mt_kahypar_free_context(context);
    mt_kahypar_free_hypergraph(hypergraph);
    mt_kahypar_free_partitioned_hypergraph(partitioned_hg);
}

static NetlistPartition read_partition_file(const Prepacker& prepacker, int num_partitions) {
    t_partition_dimension t_partition_dimensions {.x = num_partitions, .y = 1, .z = 1};
    NetlistPartition partition_map(t_partition_dimensions);
    std::ifstream partitioned_graph_file("partitioned_graph.txt");
    std::string line;

    if (!partitioned_graph_file.is_open()) {
        VPR_FATAL_ERROR(VPR_ERROR_PACK, "Could not open the partitioned graph file.\n");
    }    

    size_t count_molecules = 0;

    while (std::getline(partitioned_graph_file, line)) {
        int partition_id = std::stoi(line);
        partition_map.set_molecule_partition(PackMoleculeId(count_molecules), partition_id);
        AtomBlockId blk = prepacker.get_molecule_root_atom(PackMoleculeId(count_molecules));

        if (blk != AtomBlockId::INVALID()) {
            if (prepacker.get_expected_lowest_cost_pb_gnode(blk)->pb_type->class_type == MEMORY_CLASS) {
                partition_map.set_molecule_partition(PackMoleculeId(count_molecules), 0);
            }
        }

        count_molecules++;
    }

    return partition_map;
}

NetlistPartition NetlistPartitioner::get_graph_partitioning(int num_partitions, bool use_placement_info) {
    write_hmetis_file(use_placement_info, atom_context_, prepacker_, flat_placement_info_);

    pid_t pid = fork();
    if (pid < 0) {
        VPR_FATAL_ERROR(VPR_ERROR_PACK, "Failed to fork new process for mt-kahypar.\n");
    } else if (pid == 0) {
        call_mt_kahypar(num_partitions);
        exit(0);
    } else {
        int status;
        waitpid(pid, &status, 0);
    }

    return read_partition_file(prepacker_, num_partitions);
}

const std::vector<PackMoleculeId>& NetlistPartition::molecules(int partition) const {
    VTR_ASSERT(partition < molecules_.size());
    return molecules_[partition];
}

int NetlistPartition::get_partition(PackMoleculeId molecule) const {
    return partition_map_.at(molecule);
}

void NetlistPartition::set_molecule_partition(PackMoleculeId mol, int partition) {
    VTR_ASSERT(partition < molecules_.size());
    if(partition_map_.contains(mol)) {
        int previous_partition = partition_map_[mol];
        std::vector<PackMoleculeId>& prev_partition_vec = molecules_[previous_partition];
        prev_partition_vec.erase(std::remove(prev_partition_vec.begin(), prev_partition_vec.end(), mol), prev_partition_vec.end());
    }
    partition_map_[mol] = partition;
    molecules_[partition].push_back(mol);
}

NetlistPartition::NetlistPartition(t_partition_dimension partition_dimensions) {
    partition_dimensions_ = partition_dimensions;
    molecules_.resize(partition_dimensions.x * partition_dimensions.y * partition_dimensions.z);
}

t_partition_dimension t_partition_dimension::non_zero_shiftl(int num) {
    int new_x = this->x << num;
    int new_y = this->y << num;
    int new_z = this->z << num;

    if(new_x == 0) new_x = 1;
    if(new_y == 0) new_y = 1;
    if(new_z == 0) new_z = 1;

    return t_partition_dimension(new_x, new_y, new_z);
}

t_partition_dimension t_partition_dimension::from_index(int partition_index) {
    VTR_ASSERT(partition_index < x*y*z);
    int new_z = (partition_index % (x * y * z) ) / (x * y);
    int new_y = (partition_index % (x * y)     ) / x;
    int new_x = (partition_index % (x)         );

    return t_partition_dimension(new_x, new_y, new_z);
}

t_partition_dimension t_partition_dimension::operator/(int rhs) const {
    return t_partition_dimension(x/rhs, y/rhs, z/rhs);
}

  t_partition_dimension t_partition_dimension::operator<<(int rhs) const {
    return t_partition_dimension(x << rhs, y << rhs, z << rhs);
  }