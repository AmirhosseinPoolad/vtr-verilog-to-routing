
#include "pack.h"

#include <unordered_set>
#include "PreClusterTimingManager.h"
#include "SetupGrid.h"
#include "appack_context.h"
#include "attraction_groups.h"
#include "cluster_legalizer.h"
#include "cluster_util.h"
#include "constraints_report.h"
#include "flat_placement_types.h"
#include "globals.h"
#include "greedy_clusterer.h"
#include "partition_region.h"
#include "prepack.h"
#include "stats.h"
#include "verify_flat_placement.h"
#include "vpr_context.h"
#include "vpr_error.h"
#include "vpr_types.h"
#include "vtr_assert.h"
#include "vtr_log.h"
#include "vtr_time.h"
#include <iostream>
#include <fstream>

#include "omp.h"

#include <mtkahypar.h>
#include <flat_placement_utils.h>

#include <sys/types.h>
#include <sys/wait.h>
#include <unistd.h>
// static constexpr int thread_count = 2;

struct Position {
    float x;
    float y;
};
static void do_partitioning(int thread_count) {
    mt_kahypar_error_t error{};
    // Initialize MT-KaHyPar
    mt_kahypar_initialize(1, true);
    mt_kahypar_context_t* context = mt_kahypar_context_from_preset(DEFAULT);
    mt_kahypar_set_partitioning_parameters(context, thread_count, 0.03, CUT);
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

    std::cout << "Partitioning Results:" << std::endl;
    std::cout << "Imbalance         = " << imbalance << std::endl;
    std::cout << "Km1               = " << km1 << std::endl;
    std::cout << "Weight of Block 0 = " << block_weights[0] << std::endl;
    std::cout << "Weight of Block 1 = " << block_weights[1] << std::endl;

    mt_kahypar_free_context(context);
    mt_kahypar_free_hypergraph(hypergraph);
    mt_kahypar_free_partitioned_hypergraph(partitioned_hg);
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
                                            


static bool try_size_device_grid(const t_arch& arch,
                                 const std::map<t_logical_block_type_ptr, size_t>& num_type_instances,
                                 float target_device_utilization,
                                 const std::string& device_layout_name);

bool try_pack(const t_packer_opts& packer_opts,
              const t_analysis_opts& analysis_opts,
              const t_arch& arch,
              std::vector<t_lb_type_rr_node>* lb_type_rr_graphs,
              const Prepacker& prepacker,
              const PreClusterTimingManager& pre_cluster_timing_manager,
              const FlatPlacementInfo& flat_placement_info) {
    int thread_count = packer_opts.num_threads;
    const AtomContext& atom_ctx = g_vpr_ctx.atom();
    const DeviceContext& device_ctx = g_vpr_ctx.device();
    // The clusterer modifies the device context by increasing the size of the
    // device if needed.
    DeviceContext& mutable_device_ctx = g_vpr_ctx.mutable_device();

    std::unordered_set<AtomNetId> is_clock, is_global;
    VTR_LOG("Begin packing '%s'.\n", packer_opts.circuit_file_name.c_str());

    is_clock = alloc_and_load_is_clock();
    is_global.insert(is_clock.begin(), is_clock.end());

    size_t num_p_inputs = 0;
    size_t num_p_outputs = 0;
    for (auto blk_id : atom_ctx.netlist().blocks()) {
        auto type = atom_ctx.netlist().block_type(blk_id);
        if (type == AtomBlockType::INPAD) {
            ++num_p_inputs;
        } else if (type == AtomBlockType::OUTPAD) {
            ++num_p_outputs;
        }
    }

    VTR_LOG("\n");
    VTR_LOG("After removing unused inputs...\n");
    VTR_LOG("\ttotal blocks: %zu, total nets: %zu, total inputs: %zu, total outputs: %zu\n",
            atom_ctx.netlist().blocks().size(), atom_ctx.netlist().nets().size(), num_p_inputs, num_p_outputs);

    /* We keep attraction groups off in the first iteration,  and
     * only turn on in later iterations if some floorplan regions turn out to be overfull.
     */
    AttractionInfo attraction_groups(false);

    // We keep track of the overfilled partition regions from all pack iterations in
    // this vector. This is so that if the first iteration fails due to overfilled
    // partition regions, and it fails again, we can carry over the previous failed
    // partition regions to the current iteration.
    std::vector<PartitionRegion> overfull_partition_regions;

    t_flat_pl_loc min_coords({std::numeric_limits<float>::max(), std::numeric_limits<float>::max(), std::numeric_limits<float>::max()});
    t_flat_pl_loc max_coords({0.0f, 0.0f, 0.0f});
    t_flat_pl_loc midpoint({0.0f, 0.0f, 0.0f});

    // Verify that the Flat Placement is valid for packing.
    if (flat_placement_info.valid) {
        unsigned num_errors = verify_flat_placement_for_packing(flat_placement_info,
                                                                atom_ctx.netlist(),
                                                                prepacker);
        if (num_errors == 0) {
            VTR_LOG("Completed flat placement consistency check successfully.\n");
        } else {
            // TODO: In the future, we can just erase the flat placement and
            //       continue. It depends on what we want to happen if the
            //       flat placement is not valid.
            VPR_ERROR(VPR_ERROR_PACK,
                      "%u errors found while performing flat placement "
                      "consistency check. Aborting program.\n",
                      num_errors);
        }

        // Get min/max coordinates of the flat_placement_info
        // TODO: Find better place for this to spread out compute with locality
        for (AtomBlockId blk_id : atom_ctx.netlist().blocks()) {
            auto cur_blk_loc = flat_placement_info.get_pos(blk_id);
            max_coords.x = (cur_blk_loc.x > max_coords.x) ? cur_blk_loc.x : max_coords.x;
            max_coords.y = (cur_blk_loc.y > max_coords.y) ? cur_blk_loc.y : max_coords.y;
            max_coords.layer = (cur_blk_loc.layer > max_coords.layer) ? cur_blk_loc.layer : max_coords.layer;
            min_coords.x = (cur_blk_loc.x < min_coords.x) ? cur_blk_loc.x : min_coords.x;
            min_coords.y = (cur_blk_loc.y < min_coords.y) ? cur_blk_loc.y : min_coords.y;
            min_coords.layer = (cur_blk_loc.layer < min_coords.layer) ? cur_blk_loc.layer : min_coords.layer;
        }

        // Get the midpoint of the flat_placement_info
        midpoint = get_midpoint(min_coords, max_coords);
        VTR_LOG("Flat placement min coords: (%g, %g, %g)\n", min_coords.x, min_coords.y, min_coords.layer);
        VTR_LOG("Flat placement max coords: (%g, %g, %g)\n", max_coords.x, max_coords.y, max_coords.layer);
        VTR_LOG("Flat placement midpoint coords: (%g, %g, %g)\n", midpoint.x, midpoint.y, midpoint.layer);
        VTR_LOG("Flat placement size: (%g, %g, %g)\n", max_coords.x - min_coords.x, max_coords.y - min_coords.y, max_coords.layer - min_coords.layer);
    }

    // During clustering, a block is related to un-clustered primitives with nets.
    // This relation has three types: low fanout, high fanout, and transitive
    // high_fanout_thresholds stores the threshold for nets to a block type to
    // be considered high fanout.
    t_pack_high_fanout_thresholds high_fanout_thresholds(packer_opts.high_fanout_threshold);
    
    // Initialize the cluster legalizer.
    // Construct the APPack Context.
    APPackContext appack_ctx(flat_placement_info, device_ctx.grid);
    
    // PARALLEL TODO: add actual partitioning
    // Have options:
    // 1- Random
    // 2- Hmetis
    // 3- AP? Maybe?
    
    std::unordered_map<PackMoleculeId, int> partition_map;
    std::ofstream graph_traverse_file("graph_traverse.txt", std::ios::out);
    
    {
        vtr::ScopedFinishTimer partitioning_timer("Partitioning");
        if (!packer_opts.weighted_partitioning) {
            graph_traverse_file << g_vpr_ctx.atom().netlist().nets().size() << " " << prepacker.molecules().size() << std::endl;
        }

        // If AP is enabled, do quadrant dividing
        // Flat packer should be verified already
        if (flat_placement_info.valid && !packer_opts.weighted_partitioning) {
            // Get the number of partitions
            // Can work for thread_count of 1, 2, 4, 8, 16.
            partition_map = partition_flat_placed_mols(flat_placement_info,
                                                prepacker,
                                                thread_count,
                                                min_coords,
                                                max_coords,
                                                midpoint);
            // Print the partition map
            VTR_LOGV(packer_opts.pack_verbosity >= 3, "Partition map:\n");
            
            for (auto it = partition_map.begin(); it != partition_map.end(); ++it) {
                // auto cur_blk_loc = flat_placement_info.get_pos(it->first);
                // FIXME: Using just the first atom block id for now. Update to get entire centroid of molecule
                auto cur_blk_loc = flat_placement_info.get_pos(prepacker.get_molecule(it->first).atom_block_ids[0]);
                VTR_LOGV(packer_opts.pack_verbosity >= 3, "Molecule %zu is in partition %d ", it->first, it->second);
                VTR_LOGV(packer_opts.pack_verbosity >= 3, "With location (%g, %g, %g)\n", cur_blk_loc.x, cur_blk_loc.y, cur_blk_loc.layer);
            }
        }
        else if (true){ // Otherwise, do Hmetis graph partitioning
            std::vector<AtomBlockId> block_ids;
            std::vector<PackMoleculeId> molecule_ids;
            std::unordered_map<PackMoleculeId, std::pair<float, float>> block_positions;
            int inter_molecule_hyperedges = 0;
            float device_diameter_size = std::sqrt(device_ctx.grid.width() * device_ctx.grid.width() + device_ctx.grid.height() * device_ctx.grid.height());
            // Create graph of molecules
            bool put_space = false;
            if (graph_traverse_file.is_open()){
                for (auto net_id : g_vpr_ctx.atom().netlist().nets()){
                    molecule_ids.clear();
                    put_space = false;
                    for (auto pin_id : g_vpr_ctx.atom().netlist().net_pins(net_id)){
                        auto block_id = g_vpr_ctx.atom().netlist().pin_block(pin_id);
                        auto target_molecule_id = prepacker.get_atom_molecule(block_id);
                        auto it = std::find(molecule_ids.begin(), molecule_ids.end(), target_molecule_id);
                        if (it == molecule_ids.end()){
                            if (packer_opts.weighted_partitioning){
                                if (block_positions.find(target_molecule_id) == block_positions.end()){
                                    auto molecule_root_block_id = prepacker.get_molecule_root_atom(target_molecule_id);
                                    // Get the position of the root block in this molecule
                                    auto [block_x, block_y, block_layer] = flat_placement_info.get_pos(molecule_root_block_id);
                                    block_positions[target_molecule_id] = {block_x, block_y};
                                }
                            }
                            molecule_ids.push_back(target_molecule_id);
                        }
                    }

                    if (molecule_ids.size() > 1)
                    {
                        if (packer_opts.weighted_partitioning){
                            inter_molecule_hyperedges++;
                            auto hyperedge_weight = compute_mean_distance(block_positions);
                            graph_traverse_file << (int)(10000*(device_diameter_size  - hyperedge_weight));
                            for (auto target_molecule_id : molecule_ids){                          
                                graph_traverse_file << " ";                    
                                graph_traverse_file << ((int)target_molecule_id + 1);
                                }   
                                graph_traverse_file << "\n"; // Only write inter-molecule hyperedges
                            }                         
                        
                        else {
                            for (auto target_molecule_id : molecule_ids){
                                if (put_space){                           
                                    graph_traverse_file << " ";  
                                }                    
                                    graph_traverse_file << ((int)target_molecule_id + 1);
                                    put_space = true;
                                }       
                            } 
                    }                  
                    if (!packer_opts.weighted_partitioning){
                        graph_traverse_file << "\n";
                    }
                    block_positions.clear();       
                }
                graph_traverse_file.close();
                if (packer_opts.weighted_partitioning){
                    std::ifstream in_file("graph_traverse.txt");
                    std::ostringstream buffer;
                    buffer << in_file.rdbuf();
                    in_file.close();     
                    std::string content = buffer.str();
                    if (!content.empty() && content.back() == '\n') {
                        content.pop_back();  // remove only the last newline
                    }

                    std::ofstream graph_traverse_file("graph_traverse.txt", std::ios::out);
                    graph_traverse_file << inter_molecule_hyperedges << " " << prepacker.molecules().size() << " " << "1" << std::endl;
                    graph_traverse_file << content;
                    graph_traverse_file.close();
                }
            }
            else 
            {
                std::cerr << "Error opening file: " << std::strerror(errno) << std::endl;
            }   
            // Start partitioning
            pid_t pid = fork();
            if (pid < 0) {
                perror("fork failed");
                exit(1);
            } else if (pid == 0) {
                // In the child process: run partitioning code on its own heap
                do_partitioning(thread_count);
                exit(0);  // Ensure the child process terminates cleanly
            } else {
                // In the parent process: wait for the partitioning process to complete
                int status;
                waitpid(pid, &status, 0);
                // Now the parent can read the results (for example, reading "partitioned_graph.txt")
            }        
            // End of partitioning

            std::ifstream partitioned_graph_file("partitioned_graph.txt");
            std::string line;
            if (!partitioned_graph_file.is_open()) {
                std::cerr << "Could not open the partitioned graph file.\n";
                return 1;
            }    
            int partition_id;
            size_t count_molecules = 0;
            while (std::getline(partitioned_graph_file, line))
            {
                partition_id = std::stoi(line);
                partition_map[PackMoleculeId(count_molecules)] = partition_id;
                AtomBlockId blk = prepacker.get_molecule_root_atom(PackMoleculeId(count_molecules));
                if(blk != AtomBlockId::INVALID()) {
                    if (prepacker.get_expected_lowest_cost_pb_gnode(blk)->pb_type->class_type == MEMORY_CLASS) {
                        partition_map[PackMoleculeId(count_molecules)] = 0;
                    }
                }
                count_molecules ++;
            }
        }
        else { // Random Partitioning
            vtr::RngContainer rng(0);
            for(auto mol : prepacker.molecules()){
                partition_map[mol] = rng.irand(thread_count -1);
            }
        }

        // Re-partition all memory blocks into their own partition for AP
        // TODO: do this in a better way, this basically adds one entire O(n) stage
        if (flat_placement_info.valid) {
            for (auto blk : atom_ctx.netlist().blocks()) {
                if (prepacker.get_expected_lowest_cost_pb_gnode(blk)->pb_type->class_type == MEMORY_CLASS) {
                    partition_map[prepacker.get_atom_molecule(blk)] = 0;
                }
            }
        }

        // Verify one to one partition map
        for(auto mol : prepacker.molecules()){
            VTR_ASSERT(partition_map.contains(mol));
            VTR_ASSERT(partition_map[mol] >= 0 && partition_map[mol] <= thread_count);
        }
    }

    std::vector<std::unique_ptr<ClusterLegalizer>> cluster_legalizers;
    std::vector<std::unique_ptr<GreedyClusterer>> clusterers;
    
    for(int i = 0; i < thread_count; i++) {
        cluster_legalizers.push_back(std::make_unique<ClusterLegalizer>(atom_ctx.netlist(),
        prepacker,
        lb_type_rr_graphs,
        packer_opts.target_external_pin_util,
        high_fanout_thresholds,
        ClusterLegalizationStrategy::SKIP_INTRA_LB_ROUTE,
        packer_opts.enable_pin_feasibility_filter,
        packer_opts.pack_verbosity,
        std::ref(partition_map),
        i));

        clusterers.push_back(std::make_unique<GreedyClusterer>(packer_opts,
            analysis_opts,
            atom_ctx.netlist(),
            arch,
            high_fanout_thresholds,
            is_clock,
            is_global,
            pre_cluster_timing_manager,
            appack_ctx,
            std::ref(partition_map),
            i));
    }

    VTR_LOG("Packing with pin utilization targets: %s\n", cluster_legalizers[0]->get_target_external_pin_util().to_string().c_str());
    VTR_LOG("Packing with high fanout thresholds: %s\n", high_fanout_thresholds.to_string().c_str());


    g_vpr_ctx.mutable_atom().mutable_lookup().set_atom_pb_bimap_lock(true);
    #pragma omp parallel for num_threads(thread_count)
    for (int i = 0; i < cluster_legalizers.size(); i++) {
        double stime = omp_get_wtime();
        bool allow_unrelated_clustering = false;
        if (packer_opts.allow_unrelated_clustering == e_unrelated_clustering::ON) {
            allow_unrelated_clustering = true;
        } else if (packer_opts.allow_unrelated_clustering == e_unrelated_clustering::OFF) {
            allow_unrelated_clustering = false;
        }
    
        bool balance_block_type_util = false;
        if (packer_opts.balance_block_type_utilization == e_balance_block_type_util::ON) {
            balance_block_type_util = true;
        } else if (packer_opts.balance_block_type_utilization == e_balance_block_type_util::OFF) {
            balance_block_type_util = false;
        }
        VTR_LOG("Thread #%d\n", omp_get_thread_num());
        int pack_iteration = 1;
        while (true) {
            //Cluster the netlist
            //  num_used_type_instances: A map used to save the number of used
            //                           instances from each logical block type.
            std::map<t_logical_block_type_ptr, size_t> num_used_type_instances;
            num_used_type_instances = clusterers[i]->do_clustering(*cluster_legalizers[i],
                                                            prepacker,
                                                            allow_unrelated_clustering,
                                                            balance_block_type_util,
                                                            attraction_groups,
                                                            mutable_device_ctx);

        //Try to size/find a device
        bool fits_on_device = try_size_device_grid(arch, num_used_type_instances, packer_opts.target_device_utilization, packer_opts.device_layout);

            /* We use this bool to determine the cause for the clustering not being dense enough. If the clustering
            * is not dense enough and there are floorplan constraints, it is presumed that the constraints are the cause
            * of the floorplan not fitting, so attraction groups are turned on for later iterations.
            */
            bool floorplan_regions_overfull = floorplan_constraints_regions_overfull(overfull_partition_regions,
                                                                                    *cluster_legalizers[i],
                                                                                    device_ctx.logical_block_types);

            bool floorplan_not_fitting = (floorplan_regions_overfull || g_vpr_ctx.floorplanning().constraints.get_num_partitions() > 0);

            if (fits_on_device && !floorplan_regions_overfull) {
                break; //Done
            } else if (pack_iteration == 1 && !floorplan_not_fitting) {
                //1st pack attempt was unsuccessful (i.e. not dense enough) and we have control of unrelated clustering
                //
                //Turn it on to increase packing density
                if (packer_opts.allow_unrelated_clustering == e_unrelated_clustering::AUTO) {
                    VTR_ASSERT(allow_unrelated_clustering == false);
                    allow_unrelated_clustering = true;
                }
                if (packer_opts.balance_block_type_utilization == e_balance_block_type_util::AUTO) {
                    VTR_ASSERT(balance_block_type_util == false);
                    balance_block_type_util = true;
                }
                VTR_LOG("Packing failed to fit on device. Re-packing with: unrelated_logic_clustering=%s balance_block_type_util=%s\n",
                        (allow_unrelated_clustering ? "true" : "false"),
                        (balance_block_type_util ? "true" : "false"));
                /*
                * When running with tight floorplan constraints, some regions may become overfull with clusters (i.e.
                * the number of blocks assigned to the region exceeds the number of blocks available). When this occurs, we
                * cluster more densely to be able to adhere to the floorplan constraints. However, we do not want to cluster more
                * densely unnecessarily, as this can negatively impact wirelength. So, we have iterative approach. We check at the end
                * of every iteration if any floorplan regions are overfull. In the first iteration, we run
                * with no attraction groups (not packing more densely). If regions are overfull at the end of the first iteration,
                * we create attraction groups for partitions with overfull regions (pack those atoms more densely). We continue this way
                * until the last iteration, when we create attraction groups for every partition, if needed.
                */
            } else if (pack_iteration == 1 && floorplan_not_fitting) {
                VTR_LOG("Floorplan regions are overfull: trying to pack again using cluster attraction groups. \n");
                attraction_groups.create_att_groups_for_overfull_regions(overfull_partition_regions);
                attraction_groups.set_att_group_pulls(1);

            } else if (pack_iteration >= 2 && pack_iteration < 5 && floorplan_not_fitting) {
                if (pack_iteration == 2) {
                    VTR_LOG("Floorplan regions are overfull: trying to pack again with more attraction groups exploration. \n");
                    attraction_groups.create_att_groups_for_overfull_regions(overfull_partition_regions);
                    VTR_LOG("Pack iteration is %d\n", pack_iteration);
                } else if (pack_iteration == 3) {
                    attraction_groups.create_att_groups_for_all_regions();
                    VTR_LOG("Floorplan regions are overfull: trying to pack again with more attraction groups exploration. \n");
                    VTR_LOG("Pack iteration is %d\n", pack_iteration);
                } else if (pack_iteration == 4) {
                    attraction_groups.create_att_groups_for_all_regions();
                    VTR_LOG("Floorplan regions are overfull: trying to pack again with more attraction groups exploration and higher target pin utilization. \n");
                    VTR_LOG("Pack iteration is %d\n", pack_iteration);
                    attraction_groups.set_att_group_pulls(4);
                    t_ext_pin_util pin_util(1.0, 1.0);
                    // TODO: This line assumes the logic block name is "clb" which
                    //       may not be the case. This may need to be investigated.
                    //       Probably we should do this update of ext_pin_util for
                    //       all types that were overused. Or if that is hard, just
                    //       do it for all block types. Doing it only for a clb
                    //       string is dangerous -VB.
                    cluster_legalizers[i]->get_target_external_pin_util().set_block_pin_util("clb", pin_util);
                }

            } else { //Unable to pack densely enough: Give Up
                if (floorplan_regions_overfull) {
                    VPR_FATAL_ERROR(VPR_ERROR_OTHER,
                                    "Failed to find pack clusters densely enough to fit in the designated floorplan regions.\n"
                                    "The floorplan regions may need to be expanded to run successfully. \n");
                }

                //No suitable device found
                std::string resource_reqs;
                std::string resource_avail;
                auto& grid = g_vpr_ctx.device().grid;
                for (auto iter = num_used_type_instances.begin(); iter != num_used_type_instances.end(); ++iter) {
                    if (iter != num_used_type_instances.begin()) {
                        resource_reqs += ", ";
                        resource_avail += ", ";
                    }

                    resource_reqs += iter->first->name + ": " + std::to_string(iter->second);

                    int num_instances = 0;
                    for (auto type : iter->first->equivalent_tiles)
                        num_instances += grid.num_instances(type, -1);

                    resource_avail += iter->first->name + ": " + std::to_string(num_instances);
                }

                VPR_FATAL_ERROR(VPR_ERROR_OTHER, "Failed to find device which satisfies resource requirements required: %s (available %s)", resource_reqs.c_str(), resource_avail.c_str());
            }
        //Reset clustering for re-packing
        // PARALLEL FIXME: This probably messes up sync_netlists_to_routing_flat
        // for (auto net : g_vpr_ctx.atom().netlist().nets()) {
        //     g_vpr_ctx.mutable_atom().mutable_lookup().remove_atom_net(net);
        // }
        // g_vpr_ctx.mutable_floorplanning().cluster_constraints.clear();
        //attraction_groups.reset_attraction_groups();

            // Reset the cluster legalizer for re-clustering.
            cluster_legalizers[i]->reset();
            ++pack_iteration;
        }
        double ftime = omp_get_wtime();
        VTR_LOG("Thread %d took %f seconds in the main loop\n", i, ftime - stime);
    }

    printf("Done\n");
    /* Packing iterative improvement can be done here */
    /*       Use the re-cluster API to edit it        */
    /******************* Start *************************/
    VTR_LOG("Start the iterative improvement process\n");
    //iteratively_improve_packing(*packer_opts, clustering_data, 2);
    VTR_LOG("the iterative improvement process is done\n");

    /*
     * auto& cluster_ctx = g_vpr_ctx.clustering();
     * for (auto& blk_id : g_vpr_ctx.clustering().clb_nlist.blocks()) {
     * free_pb_stats_recursive(cluster_ctx.clb_nlist.block_pb(blk_id));
     * }
     */
    /******************** End **************************/
    g_vpr_ctx.mutable_atom().mutable_lookup().set_atom_pb_bimap_lock(false);

    verify_clustering(cluster_legalizers);
    AtomPBBimap final_atom_pb;
    for (auto atom_blk : g_vpr_ctx.atom().netlist().blocks()){
        for(auto &cluster_legalizer : cluster_legalizers) {
            auto atom_pb = cluster_legalizer->atom_pb_lookup().atom_pb(atom_blk);
            if (atom_pb != nullptr) {
                final_atom_pb.set_atom_pb(atom_blk, atom_pb);
                break;
            }
        }
    }
    
    g_vpr_ctx.mutable_atom().mutable_lookup().set_atom_to_pb_bimap(final_atom_pb);

    for (auto atom_blk : g_vpr_ctx.atom().netlist().blocks()){
        VTR_ASSERT(g_vpr_ctx.mutable_atom().mutable_lookup().atom_pb_bimap().atom_pb(atom_blk) != nullptr);
    }

    //check clustering and output it
    // cluster_legalizers[0]->mutable_atom_pb_lookup() = final_atom_pb_bimap;
    check_and_output_clustering(cluster_legalizers, packer_opts, is_clock, &arch);
    VTR_LOG("\n");
    VTR_LOG("Netlist conversion complete.\n");
    VTR_LOG("\n");

    return true;
}

std::unordered_set<AtomNetId> alloc_and_load_is_clock() {
    /* Looks through all the atom blocks to find and mark all the clocks, by setting
     * the corresponding entry by adding the clock to is_clock.
     * only for an error check.                                                */

    std::unordered_set<AtomNetId> is_clock;

    /* Want to identify all the clock nets.  */
    auto& atom_ctx = g_vpr_ctx.atom();

    for (auto blk_id : atom_ctx.netlist().blocks()) {
        for (auto pin_id : atom_ctx.netlist().block_clock_pins(blk_id)) {
            auto net_id = atom_ctx.netlist().pin_net(pin_id);
            if (!is_clock.count(net_id)) {
                is_clock.insert(net_id);
            }
        }
    }

    return (is_clock);
}

static bool try_size_device_grid(const t_arch& arch,
                                 const std::map<t_logical_block_type_ptr, size_t>& num_type_instances,
                                 float target_device_utilization,
                                 const std::string& device_layout_name) {
    auto& device_ctx = g_vpr_ctx.mutable_device();

    //Build the device
    auto grid = create_device_grid(device_layout_name, arch.grid_layouts, num_type_instances, target_device_utilization);

    /*
     *Report on the device
     */
    VTR_LOG("FPGA sized to %zu x %zu (%s)\n", grid.width(), grid.height(), grid.name().c_str());

    bool fits_on_device = true;

    float device_utilization = calculate_device_utilization(grid, num_type_instances);
    VTR_LOG("Device Utilization: %.2f (target %.2f)\n", device_utilization, target_device_utilization);
    std::map<t_logical_block_type_ptr, float> type_util;
    for (const auto& type : device_ctx.logical_block_types) {
        if (is_empty_type(&type)) continue;

        auto itr = num_type_instances.find(&type);
        if (itr == num_type_instances.end()) continue;

        float num_instances = itr->second;
        float util = 0.;

        float num_total_instances = 0.;
        for (const auto& equivalent_tile : type.equivalent_tiles) {
            num_total_instances += device_ctx.grid.num_instances(equivalent_tile, -1);
        }

        if (num_total_instances != 0) {
            util = num_instances / num_total_instances;
        }
        type_util[&type] = util;

        if (util > 1.) {
            fits_on_device = false;
        }
        VTR_LOG("\tBlock Utilization: %.2f Type: %s\n", util, type.name.c_str());
    }
    VTR_LOG("\n");

    return fits_on_device;
}

std::unordered_map<PackMoleculeId, int> partition_flat_placed_mols(const FlatPlacementInfo& flat_placement_info,
    const Prepacker& prepacker,
    int thread_count,
    const t_flat_pl_loc& min_coords,
    const t_flat_pl_loc& max_coords,
    const t_flat_pl_loc& midpoint) {
    /* This function partitions the region into smaller rectangles and assigns each molecule to a partition.
     * The partition is determined by the location of the molecule in the region. */

    std::unordered_map<PackMoleculeId, int> partition_map;
    
    // Get the number of partitions
    // Can work for thread_count of 1, 2, 4, 8, 16.
    int num_partitions = thread_count;
    int num_partitions_x, num_partitions_y, num_partitions_layer = 0;//sqrt(num_partitions);

    if (num_partitions == 1) {
        // No need to partition, just return 0
        for (auto mol : prepacker.molecules()) {
            partition_map[mol] = 0;
        }
        return partition_map; 
    }

    // Hardcode thread count to partitions for now
    switch(num_partitions) {
        case 2:
            num_partitions_x = 2;
            num_partitions_y = 1;
            break;
        case 4:
            num_partitions_x = 2;
            num_partitions_y = 2;
            break;
        case 8:
            num_partitions_x = 4;
            num_partitions_y = 2;
            break;
        case 16:
            num_partitions_x = 4;
            num_partitions_y = 4;
            break;
        default:
            VPR_FATAL_ERROR(VPR_ERROR_PACK, "Unsupported number of partitions: %d\n", num_partitions);
    }


    // Not that many layers, should be able to just use the number of layers as the partition size
    // TODO: For now, number of layer paritions is hardcoded to 1
    num_partitions_layer = 1;
    int partition_size_layer = 1; //(max_coords.layer - min_coords.layer) / num_partitions_layer; 
    // Increase max coords by 1% to make max_coord block be num_partitions_xy-1.
    double partition_size_x = (max_coords.x - min_coords.x) / num_partitions_x;
    double partition_size_y = (max_coords.y - min_coords.y) / num_partitions_y;
    VTR_LOG("Partition size: (%lf, %lf, %d)\n", partition_size_x, partition_size_y, partition_size_layer);

    // Combine partitions until the total number of partitions is equal to the number of threads
    // This is done by combining the x and y partitions
    // The number of partitions is the product of the number of partitions in each dimension
    // for 


    // Get the partition for each molecule
    for (auto mol : prepacker.molecules()) {
        // auto blk_id = flat_placement_info.blx_id[mol];
        // auto cur_blk_loc = flat_placement_info.get_pos(mol);
        auto cur_blk_loc = flat_placement_info.get_pos(prepacker.get_molecule(mol).atom_block_ids[0]);

        auto x = cur_blk_loc.x;
        auto y = cur_blk_loc.y;
        auto layer = cur_blk_loc.layer;

        // Get the partition for the molecule
        int partition_x = (int)((x - min_coords.x) / partition_size_x);
        int partition_y = (int)((y - min_coords.y) / partition_size_y);
        partition_x = partition_x == num_partitions_x ? partition_x - 1 : partition_x; //
        partition_y = partition_y == num_partitions_y ? partition_y - 1 : partition_y;

        // // Move blocks close to the right border of partition to the next partiton
        // // TODO: only works when mem/dsp/etc are arranged in columns. find a better solution
        // constexpr int border_threshold = 1;
        // if (partition_x != num_partitions_x - 1) {
        //     if (std::abs(x - (partition_size_x * (partition_x + 1) + min_coords.x)) <= border_threshold) {
        //         partition_x++;
        //     }
        // }
        // TODO:Always 1 layer for now
        int partition_layer = 0;// (int)((layer - min_coords.layer) / partition_size_layer);

        // Assign the molecule to the partition
        int partition_id = partition_x + (partition_y * num_partitions_x); // + partition_layer * num_partitions_x * num_partitions_y;

        // Check if the partition is valid. This should never happen
        if (partition_x < 0 || partition_x >= num_partitions_x ||
            partition_y < 0 || partition_y >= num_partitions_y ||
            partition_layer < 0 || partition_layer >= num_partitions_layer ||
            partition_id < 0 || partition_id >= num_partitions) {
            VPR_FATAL_ERROR(VPR_ERROR_PACK, "Partition: (%d, %d, %d) = (%d) is out of bounds: max partitions(%d, %d, %d) = (%d)\n", partition_x, partition_y, partition_layer, partition_id, num_partitions_x, num_partitions_y, num_partitions_layer, num_partitions);
        }
        //VTR_LOG("(%f, %f) = (%d)\n", x, y, partition_id);
        partition_map[mol] = partition_id;
    }
    return partition_map;
}