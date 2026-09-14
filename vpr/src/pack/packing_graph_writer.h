#pragma once

#include <string>
#include <unordered_set>
#include <vector>

#include "atom_netlist_fwd.h"
#include "lb_type_rr_node_types.h"

struct t_intra_lb_net;
struct t_mode_selection_status;

/**
 * @brief Append the router's current atom-level cluster graph and terminal
 *        routing result to packing_graph.json.
 */
void write_packing_graph(const std::string& cluster_type,
                         const std::unordered_set<AtomBlockId>& atoms,
                         const std::vector<t_intra_lb_net>& intra_lb_nets,
                         const std::vector<t_lb_type_rr_node>& lb_graph,
                         bool is_impossible,
                         bool route_succeeded,
                         const t_mode_selection_status& mode_status);
