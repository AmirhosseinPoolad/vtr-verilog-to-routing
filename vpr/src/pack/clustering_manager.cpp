#include "clustering_manager.h"
#include "atom_netlist.h"
#include "cluster_legalizer.h"
#include "vpr_types.h"


/*
 * @brief Check the atom blocks of a cluster pb. Used in the verify method.
 */
/* TODO: May want to check that all atom blocks are actually reached */
static void check_cluster_atom_blocks(t_pb* pb, std::unordered_set<AtomBlockId>& blocks_checked, const AtomPBBimap& atom_pb_lookup, const AtomNetlist& atom_netlist) {
    const t_pb_type* pb_type = pb->pb_graph_node->pb_type;
    if (pb_type->is_primitive()) {
        /* primitive */
        AtomBlockId blk_id = atom_pb_lookup.pb_atom(pb);
        if (blk_id) {
            if (blocks_checked.count(blk_id)) {
                VPR_FATAL_ERROR(VPR_ERROR_PACK,
                                "pb %s contains atom block %s but atom block is already contained in another pb.\n",
                                pb->name, atom_netlist.block_name(blk_id).c_str());
            }
            blocks_checked.insert(blk_id);
            if (pb != atom_pb_lookup.atom_pb(blk_id)) {
                VPR_FATAL_ERROR(VPR_ERROR_PACK,
                                "pb %s contains atom block %s but atom block does not link to pb.\n",
                                pb->name, atom_netlist.block_name(blk_id).c_str());
            }
        }
    } else {
        /* this is a container pb, all container pbs must contain children */
        bool has_child = false;
        for (int i = 0; i < pb_type->modes[pb->mode].num_pb_type_children; i++) {
            for (int j = 0; j < pb_type->modes[pb->mode].pb_type_children[i].num_pb; j++) {
                if (pb->child_pbs[i] != nullptr) {
                    if (pb->child_pbs[i][j].name != nullptr) {
                        has_child = true;
                        check_cluster_atom_blocks(&pb->child_pbs[i][j], blocks_checked, atom_pb_lookup, atom_netlist);
                    }
                }
            }
        }
        VTR_ASSERT(has_child);
    }
}

void ClusteringManager::verify_clustering() const {
    std::unordered_set<AtomBlockId> atoms_checked;
    int clustered_mols = 0;
    for (auto& cluster_legalizer_ptr : cluster_legalizers_) {
        clustered_mols += cluster_legalizer_ptr->clusters().size();
    }
    if (clustered_mols == 0) {
        VTR_LOG_WARN("Packing produced no clustered blocks");
    }
    /*
    * Check that each atom block connects to one physical primitive and that the primitive links up to the parent clb
    */
    for (auto blk_id : atom_netlist_.blocks()) {
        //Each atom should be part of a pb
        const t_pb* atom_pb = nullptr;

        // Find the ClusterLegalizer object that has packed blk_id
        size_t cluster_legalizer_index = -1;
        for(size_t i = 0; i < cluster_legalizers_.size(); i++) {
            const t_pb* pb = cluster_legalizers_[i]->atom_pb_lookup().atom_pb(blk_id);
            if(pb != nullptr) {
                VTR_ASSERT_MSG(atom_pb == nullptr, "Atom is packed in two different ClusterLegalizers");
                atom_pb = pb;
                cluster_legalizer_index = i;
            }
        }
        const ClusterLegalizer& cluster_legalizer = *cluster_legalizers_[cluster_legalizer_index];

        if (!atom_pb) {
            VPR_FATAL_ERROR(VPR_ERROR_PACK,
                            "Atom block %s with ID %d is not mapped to a pb\n",
                            atom_netlist_.block_name(blk_id).c_str(), (int)blk_id);
        }

        //Check the reverse mapping is consistent
        if (cluster_legalizer.atom_pb_lookup().pb_atom(atom_pb) != blk_id) {
            VPR_FATAL_ERROR(VPR_ERROR_PACK,
                            "pb %s does not contain atom block %s but atom block %s maps to pb.\n",
                            atom_pb->name,
                            atom_netlist_.block_name(blk_id).c_str(),
                            atom_netlist_.block_name(blk_id).c_str());
        }

        VTR_ASSERT(atom_netlist_.block_name(blk_id) == atom_pb->name);

        const t_pb* cur_pb = atom_pb;
        while (cur_pb->parent_pb) {
            cur_pb = cur_pb->parent_pb;
            VTR_ASSERT(cur_pb->name);
        }

        LegalizationClusterId cluster_id = cluster_legalizer.get_atom_cluster(blk_id);
        if (cluster_id == LegalizationClusterId::INVALID()) {
            VPR_FATAL_ERROR(VPR_ERROR_PACK,
                            "Atom %s is not mapped to a CLB\n",
                            atom_netlist_.block_name(blk_id).c_str());
        }

        if (cur_pb != cluster_legalizer.get_cluster_pb(cluster_id)) {
            VPR_FATAL_ERROR(VPR_ERROR_PACK,
                            "CLB %s does not match CLB contained by pb %s.\n",
                            cur_pb->name, atom_pb->name);
        }
    }

    /* Check that I do not have spurious links in children pbs */
    for(size_t i = 0; i < cluster_legalizers_.size(); i++) {
        for (LegalizationClusterId cluster_id : cluster_legalizers_[i]->clusters()) {
            if (!cluster_id.is_valid())
                continue;
            check_cluster_atom_blocks(cluster_legalizers_[i]->get_cluster_pb(cluster_id), atoms_checked, cluster_legalizers_[i]->atom_pb_lookup(), atom_netlist_);
        }
    }

    for (auto blk_id : atom_netlist_.blocks()) {
        if (!atoms_checked.contains(blk_id)) {
            VPR_FATAL_ERROR(VPR_ERROR_PACK,
                            "Atom block %s not found in any cluster.\n",
                            atom_netlist_.block_name(blk_id).c_str());
        }
    }
}

ClusteringManager::ClusteringManager(const t_packer_opts& packer_opts,
                                     const t_analysis_opts& analysis_opts,
                                     const AtomNetlist& atom_netlist,
                                     const t_arch& arch,
                                     const Prepacker& prepacker,
                                     std::vector<t_lb_type_rr_node>* lb_type_rr_graphs,
                                     const t_pack_high_fanout_thresholds& high_fanout_thresholds,
                                     ClusterLegalizationStrategy cluster_legalization_strategy,
                                     bool enable_pin_feasibility_filter,
                                     const LogicalModels& models,
                                     const std::unordered_set<AtomNetId>& is_clock,
                                     const std::unordered_set<AtomNetId>& is_global,
                                     const PreClusterTimingManager& pre_cluster_timing_manager,
                                     const APPackContext& appack_ctx,
                                     int log_verbosity,
                                     const NetlistPartition& netlist_partition) :
                                     atom_netlist_(atom_netlist) {
    int thread_count = packer_opts.num_threads;
    target_external_pin_util_ = t_ext_pin_util_targets(packer_opts.target_external_pin_util);
    t_pack_high_fanout_thresholds high_fanout_thresholds_(packer_opts.high_fanout_threshold);

    for (int i = 0; i < thread_count; i++) {
        cluster_legalizers_.push_back(std::make_unique<ClusterLegalizer>(atom_netlist,
        prepacker,
        lb_type_rr_graphs,
        target_external_pin_util_,
        high_fanout_thresholds,
        cluster_legalization_strategy,
        enable_pin_feasibility_filter,
        models,
        log_verbosity,
        netlist_partition));

        clusterers_.push_back(std::make_unique<GreedyClusterer>(packer_opts,
            analysis_opts,
            atom_netlist,
            arch,
            high_fanout_thresholds,
            is_clock,
            is_global,
            pre_cluster_timing_manager,
            appack_ctx,
            netlist_partition));
    }
}