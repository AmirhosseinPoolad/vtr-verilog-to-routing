#pragma once

#include <memory>
#include <optional>
#include <vector>
#include "atom_netlist.h"
#include "atom_netlist_fwd.h"
#include "atom_pb_bimap.h"
#include "cluster_legalizer.h"
#include "greedy_clusterer.h"
#include "vpr_types.h"
#include "vtr_assert.h"

class ClusteringManager {

// iterators that might be needed:
// cluster_legalizers.clusters

public:
    /// @brief Returns a constant reference to the target_external_pin_util object.
    inline const t_ext_pin_util_targets& get_target_external_pin_util() const {
        return target_external_pin_util_;
    }

    inline t_ext_pin_util_targets& get_mut_target_external_pin_util() {
        return target_external_pin_util_;
    }

    const std::vector<std::unique_ptr<ClusterLegalizer>>& cluster_legalizers() const {
        return cluster_legalizers_;
    }

    const std::vector<std::unique_ptr<GreedyClusterer>>& clusterers() const {
        return clusterers_;
    }

    size_t clustering_classes_count(){
        VTR_ASSERT_SAFE(cluster_legalizers_.size() == clusterers_.size());
        return cluster_legalizers_.size();
    }

    AtomPBBimap atom_to_pb(){
        AtomPBBimap final_atom_pb;
        for (auto atom_blk : atom_netlist_.blocks()) {
            for (auto& cluster_legalizer : cluster_legalizers_) {
                auto atom_pb = cluster_legalizer->atom_pb_lookup().atom_pb(atom_blk);
                if (atom_pb != nullptr) {
                    final_atom_pb.set_atom_pb(atom_blk, atom_pb);
                    break;
                }
            }
        }
        return final_atom_pb;
    }

    void reset_clustering_state() {
        for(auto& cluster_legalizer : cluster_legalizers_) {
            cluster_legalizer->reset();
        }
    }

    std::pair<const ClusterLegalizer&, bool> get_atom_cluster_legalizer(AtomBlockId blk_id) const {
        for(const std::unique_ptr<ClusterLegalizer>& cluster_legalizer_ptr : cluster_legalizers_) {
            if (cluster_legalizer_ptr->is_atom_clustered(blk_id)) {
                return {*cluster_legalizer_ptr, true};
            }
        }
        return {*cluster_legalizers_[0], false};
    }

    /*
     * @brief Verify that all atoms have been clustered into some cluster.
     *
     * This will not verify if all the clusters are fully legal.
     */
    void verify_clustering() const;

    ClusteringManager(const t_packer_opts& packer_opts,
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
                                     const NetlistPartition& netlist_partition);

    ClusteringManager();

private:
    std::vector<std::unique_ptr<ClusterLegalizer>> cluster_legalizers_;
    std::vector<std::unique_ptr<GreedyClusterer>> clusterers_;

    /// @brief The maximum fractional utilization of cluster external
    ///        input/output pins during packing (between 0 and 1).
    t_ext_pin_util_targets target_external_pin_util_;

    t_pack_high_fanout_thresholds high_fanout_thresholds_;//(packer_opts.high_fanout_threshold);

    const AtomNetlist& atom_netlist_;
};