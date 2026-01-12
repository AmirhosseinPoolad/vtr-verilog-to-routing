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
#include "vtr_assert.h"
#include "vtr_log.h"
#include "vtr_random.h"
#include "router_lookahead_map_utils.h"

static constexpr int num_landmarks = 112;

DifferentialLookahead::DifferentialLookahead(const t_det_routing_arch& det_routing_arch, bool is_flat, int route_verbosity) {
    map_lookahead = std::make_unique<MapLookahead>(det_routing_arch, is_flat, route_verbosity);
    const auto& device = g_vpr_ctx.device();
    const auto& rr_graph = device.rr_graph;
    const auto& spatial_lookup = device.rr_graph.node_lookup();
    // select landmark nodes
    std::vector<std::pair<int, int>> landmark_locs;
    vtr::RandomNumberGenerator rng;

    auto [x_size, y_size, z_size] = device.grid.dim_sizes();
    // num_landmarks / 4 in each side

    // bottom
    for (int i = 0; i < num_landmarks / 4; i++) {
        landmark_locs.push_back({i * (x_size / (num_landmarks / 4)), 0});
    }

    // top
    for (int i = 0; i < num_landmarks / 4; i++) {
        landmark_locs.push_back({i * (x_size / (num_landmarks / 4)), y_size - 1});
    }

    // left
    for (int i = 0; i < num_landmarks / 4; i++) {
        landmark_locs.push_back({0, i * (y_size / (num_landmarks / 4))});
    }

    // right
    for (int i = 0; i < num_landmarks / 4; i++) {
        landmark_locs.push_back({x_size, i * (y_size / (num_landmarks / 4))});
    }

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
    float delay_cost = std::numeric_limits<float>::lowest();

    for (int i = 0; i < num_landmarks; i++) {
        delay_cost = std::max(delay_cost, std::abs(landmark_costs[i][(size_t)node] - landmark_costs[i][(size_t)target_node]));
    }
    auto map_res = map_lookahead->get_expected_delay_and_cong(node, target_node, params, R_upstream);
    return {delay_cost * params.criticality, map_res.second};
}
