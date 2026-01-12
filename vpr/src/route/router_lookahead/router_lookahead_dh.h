#pragma once

#include <memory>
#include <vector>
#include "router_lookahead_map.h"
#include "vpr_types.h"
#include "vpr_error.h"
#include "router_lookahead.h"
#include "router_lookahead_map_utils.h"

class DifferentialLookahead : public RouterLookahead {
  public:
    explicit DifferentialLookahead(const t_det_routing_arch& det_routing_arch, bool is_flat, int route_verbosity);

  private:
  std::vector<std::vector<float>> landmark_costs;
  std::unique_ptr<MapLookahead> map_lookahead;

  protected:
    float get_expected_cost(RRNodeId node, RRNodeId target_node, const t_conn_cost_params& params, float R_upstream) const override;
    std::pair<float, float> get_expected_delay_and_cong(RRNodeId node, RRNodeId target_node, const t_conn_cost_params& params, float R_upstream) const override;

    void compute(const std::vector<t_segment_inf>& segment_inf) override {
        map_lookahead->compute(segment_inf);
    }

    void compute_intra_tile() override {
        VPR_THROW(VPR_ERROR_ROUTE, "ClassicLookahead::compute_intra_time unimplemented");
    }

    void read(const std::string& /*file*/) override {
        VPR_THROW(VPR_ERROR_ROUTE, "Read not supported for NoOpLookahead");
    }

    void read_intra_cluster(const std::string& /*file*/) override {
        VPR_THROW(VPR_ERROR_ROUTE, "read_intra_cluster not supported for NoOpLookahead");
    }

    void write(const std::string& /*file*/) const override {
        VPR_THROW(VPR_ERROR_ROUTE, "Write not supported for NoOpLookahead");
    }

    void write_intra_cluster(const std::string& /*file*/) const override {
        VPR_THROW(VPR_ERROR_ROUTE, "write_intra_cluster not supported for NoOpLookahead");
    }

    float get_opin_distance_min_delay(int /*physical_tile_idx*/, int /*from_layer*/, int /*to_layer*/, int /*dx*/, int /*dy*/) const override {
        return -1.;
    }
};