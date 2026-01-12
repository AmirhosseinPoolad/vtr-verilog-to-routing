#include "router_lookahead_dh.h"
#include "connection_router_interface.h"
#include <cstddef>
#include <cstdlib>
#include <limits>
#include <memory>
#include <vector>
#include "globals.h"
#include "router_lookahead_map.h"
#include "rr_graph_fwd.h"
#include "rr_node_types.h"
#include "vtr_assert.h"
#include "vtr_log.h"
#include "vtr_random.h"
#include "router_lookahead_map_utils.h"

static constexpr int num_landmarks = 112;

static std::vector<std::pair<int, int>> get_uniform_perimeter_points(int k, int x_size, int y_size) {
    std::vector<std::pair<int, int>> points;
    
    // Safety check: ensure the inner rectangle exists
    if (k <= 0 || x_size <= 2 || y_size <= 2) {
        return points;
    }

    // Coordinates of the inner rectangle boundaries
    int x_min = 1, y_min = 1;
    int x_max = x_size - 2;
    int y_max = y_size - 2;

    int w = x_max - x_min;
    int h = y_max - y_min;
    long long perimeter = 2LL * (w + h);

    for (int i = 0; i < k; ++i) {
        // Distribute k points along the perimeter using integer scaling
        long long d = (i * perimeter) / k;
        int x, y;

        if (d < w) {
            // Bottom edge: moving right
            x = x_min + (int)d;
            y = y_min;
        } else if (d < (w + h)) {
            // Right edge: moving up
            x = x_max;
            y = y_min + (int)(d - w);
        } else if (d < (2LL * w + h)) {
            // Top edge: moving left
            x = x_max - (int)(d - (w + h));
            y = y_max;
        } else {
            // Left edge: moving down
            x = x_min;
            y = y_max - (int)(d - (2LL * w + h));
        }
        
        points.push_back({x, y});
    }

    return points;
}

DifferentialLookahead::DifferentialLookahead(const t_det_routing_arch& det_routing_arch, bool is_flat, int route_verbosity) {
    map_lookahead = std::make_unique<MapLookahead>(det_routing_arch, is_flat, route_verbosity);
    const auto& device = g_vpr_ctx.device();
    const auto& rr_graph = device.rr_graph;
    const auto& spatial_lookup = device.rr_graph.node_lookup();
    // select landmark nodes
    std::vector<std::pair<int, int>> landmark_locs;
    vtr::RandomNumberGenerator rng;

    auto [z_size, x_size, y_size] = device.grid.dim_sizes();
    // num_landmarks / 4 in each side

    landmark_locs = get_uniform_perimeter_points(num_landmarks, x_size, y_size);
    VTR_LOG("%d %d %d\n", landmark_locs.size(), x_size, y_size);
    VTR_ASSERT(landmark_locs.size() == num_landmarks);

    // initialize the size of the landmark costs vectors
    size_t num_nodes = rr_graph.num_nodes();

    landmark_costs.resize(num_landmarks);

    for (auto& landmark_cost : landmark_costs) {
        landmark_cost.resize(num_nodes);
    }

    VTR_LOG("Landmark costs memory cost (bytes): %d\n", num_landmarks * num_nodes * sizeof(float));

    // start single-source all-sinks routings from all landmark locations
    for (int landmark_index = 0; landmark_index < num_landmarks; landmark_index++) {
        vtr::vector<RRNodeId, bool> node_expanded;
        node_expanded.resize(rr_graph.num_nodes());
        std::fill(node_expanded.begin(), node_expanded.end(), false);

        vtr::vector<RRNodeId, float> node_visited_costs;
        node_visited_costs.resize(rr_graph.num_nodes());
        std::fill(node_visited_costs.begin(), node_visited_costs.end(), -1.0);

        const auto landmark_loc = landmark_locs[landmark_index];

        // A priority queue for expansion
        std::priority_queue<util::PQ_Entry> pq;

        std::vector<RRNodeId> starting_nodes;
        // add these to the pq
        auto nodes = spatial_lookup.find_channel_nodes(0, landmark_loc.first, landmark_loc.second, e_rr_type::CHANX);
        starting_nodes.insert(starting_nodes.end(), nodes.begin(), nodes.end());

        nodes =spatial_lookup.find_channel_nodes(0, landmark_loc.first, landmark_loc.second, e_rr_type::CHANY);
        starting_nodes.insert(starting_nodes.end(), nodes.begin(), nodes.end());

        nodes =spatial_lookup.find_channel_nodes(0, landmark_loc.first + 1, landmark_loc.second, e_rr_type::CHANX);
        starting_nodes.insert(starting_nodes.end(), nodes.begin(), nodes.end());

        nodes =spatial_lookup.find_channel_nodes(0, landmark_loc.first, landmark_loc.second + 1, e_rr_type::CHANY);
        starting_nodes.insert(starting_nodes.end(), nodes.begin(), nodes.end());

        // First entry has no upstream delay or congestion
        for (auto node : starting_nodes) {
            pq.emplace(node, UNDEFINED, 0, 0, 0, true);
        }
        // Now do routing
        while (!pq.empty()) {
            util::PQ_Entry current = pq.top();
            pq.pop();

            RRNodeId curr_node = current.rr_node;

            // Check that we haven't already expanded from this node
            if (node_expanded[curr_node]) {
                continue;
            }

            landmark_costs[landmark_index][(size_t)curr_node] = current.cost;

            // expand node
            for (t_edge_size edge : rr_graph.edges(curr_node)) {
                RRNodeId child_node = rr_graph.edge_sink_node(curr_node, edge);
                int switch_ind = size_t(rr_graph.edge_switch(curr_node, edge));

                if (node_expanded[child_node]) {
                    continue;
                }

                util::PQ_Entry child_entry(child_node, switch_ind, current.delay,
                                           current.R_upstream, current.congestion_upstream, false);

                //VTR_ASSERT(child_entry.cost >= 0); //Assertion fails in practise. TODO: debug

                /* skip this child if it has been visited with smaller cost */
                if (node_visited_costs[child_node] >= 0 && node_visited_costs[child_node] < child_entry.cost) {
                    continue;
                }

                /* finally, record the cost with which the child was visited and put the child entry on the queue */
                node_visited_costs[child_node] = child_entry.cost;
                pq.push(child_entry);
            }

            node_expanded[curr_node] = true;
        }
    }
}

float DifferentialLookahead::get_expected_cost(RRNodeId node, RRNodeId target_node, const t_conn_cost_params& params, float R_upstream) const {
    auto costs = get_expected_delay_and_cong(node, target_node, params, R_upstream);
    return costs.first + costs.second;
}
std::pair<float, float> DifferentialLookahead::get_expected_delay_and_cong(RRNodeId node, RRNodeId target_node, const t_conn_cost_params& params, float R_upstream) const {
    auto map_res = map_lookahead->get_expected_delay_and_cong(node, target_node, params, R_upstream);

    const auto& device_ctx = g_vpr_ctx.device();
    const auto& rr_graph = device_ctx.rr_graph;
    if (rr_graph.node_type(node) != e_rr_type::CHANZ && rr_graph.node_type(node) != e_rr_type::CHANX && rr_graph.node_type(node) != e_rr_type::CHANY) {
        return map_res;
    }

    float delay_cost = std::numeric_limits<float>::lowest();

    for (int i = 0; i < num_landmarks; i++) {
        delay_cost = std::max(delay_cost, std::abs(landmark_costs[i][(size_t)node] - landmark_costs[i][(size_t)target_node]));
    }

    return {delay_cost * params.criticality, map_res.second};
}
