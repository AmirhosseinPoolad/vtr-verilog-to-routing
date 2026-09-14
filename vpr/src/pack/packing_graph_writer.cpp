#include "packing_graph_writer.h"

#include <algorithm>
#include <fstream>
#include <iomanip>
#include <sstream>
#include <string>
#include <string_view>
#include <unordered_map>
#include <unordered_set>
#include <vector>

#include "atom_netlist.h"
#include "cluster_router.h"
#include "globals.h"
#include "physical_types.h"
#include "vpr_context.h"
#include "vtr_assert.h"
#include "vtr_log.h"

namespace {

constexpr size_t EXTERNAL_INPUT_NODE = 0;
constexpr size_t EXTERNAL_OUTPUT_NODE = 1;
constexpr size_t FIRST_ATOM_NODE = 2;

/**
 * @brief Atom-node attributes serialized for an intra-cluster packing graph.
 */
struct t_packing_graph_atom_node {
    size_t atom_id;    ///< Atom netlist block ID.
    std::string model; ///< Atom logical model name.
    int primitive_num; ///< Cluster-relative physical primitive number.
};

/**
 * @brief Directed edge in an intra-cluster packing graph.
 */
struct t_packing_graph_edge {
    size_t source_node; ///< Dense source-node ID.
    size_t sink_node;   ///< Dense sink-node ID.
    size_t atom_net_id; ///< Atom netlist net ID associated with the edge.
};

/**
 * @brief One labeled intra-cluster packing graph sample.
 */
struct t_packing_graph_record {
    std::string cluster_type;                          ///< Logical block type being routed.
    bool is_impossible;                                ///< A required connection had no legal path.
    bool route_succeeded;                              ///< The final routing iteration succeeded.
    std::vector<t_packing_graph_atom_node> atom_nodes; ///< Clustered atom nodes.
    std::vector<t_packing_graph_edge> edges;           ///< Expanded directed atom connections.
};

/**
 * @brief Escape and quote a string for JSON output.
 */
std::string json_string(const std::string_view value) {
    std::ostringstream output;
    output << '"';
    for (unsigned char character : value) {
        switch (character) {
            case '"':
                output << "\\\"";
                break;
            case '\\':
                output << "\\\\";
                break;
            case '\b':
                output << "\\b";
                break;
            case '\f':
                output << "\\f";
                break;
            case '\n':
                output << "\\n";
                break;
            case '\r':
                output << "\\r";
                break;
            case '\t':
                output << "\\t";
                break;
            default:
                if (character < 0x20) {
                    output << "\\u"
                           << std::hex << std::setw(4) << std::setfill('0')
                           << static_cast<int>(character) << std::dec;
                } else {
                    output << static_cast<char>(character);
                }
        }
    }
    output << '"';
    return output.str();
}

/**
 * @brief Find the primitive number associated with each terminal RR node.
 */
std::unordered_map<int, int> terminal_primitive_nums(const std::vector<t_lb_type_rr_node>& lb_graph) {
    std::unordered_map<int, int> primitive_nums;

    for (size_t inode = 0; inode < lb_graph.size(); ++inode) {
        const t_lb_type_rr_node& rr_node = lb_graph[inode];
        if (rr_node.pb_graph_pin == nullptr
            || !rr_node.pb_graph_pin->parent_node->pb_type->is_primitive()) {
            continue;
        }

        int primitive_num = rr_node.pb_graph_pin->parent_node->primitive_num;
        primitive_nums.insert({static_cast<int>(inode), primitive_num});

        for (int imode = 0; imode < rr_node.num_modes; ++imode) {
            for (int iedge = 0; iedge < rr_node.num_fanout[imode]; ++iedge) {
                int sink_node = rr_node.outedges[imode][iedge].node_index;
                if (lb_graph[sink_node].type == LB_SINK) {
                    primitive_nums.insert({sink_node, primitive_num});
                }
            }
        }
    }

    return primitive_nums;
}

/**
 * @brief Build a serializable graph record from the router's current state.
 */
t_packing_graph_record make_packing_graph(const std::string& cluster_type,
                                          const std::unordered_set<AtomBlockId>& atoms_added,
                                          const std::vector<t_intra_lb_net>& intra_lb_nets,
                                          const std::vector<t_lb_type_rr_node>& lb_graph,
                                          bool is_impossible,
                                          bool route_succeeded) {
    const AtomContext& atom_ctx = g_vpr_ctx.atom();
    const AtomNetlist& atom_netlist = atom_ctx.netlist();
    const LogicalModels& models = g_vpr_ctx.device().arch->models;

    t_packing_graph_record graph;
    graph.cluster_type = cluster_type;
    graph.is_impossible = is_impossible;
    graph.route_succeeded = route_succeeded;

    std::vector<AtomBlockId> atoms(atoms_added.begin(), atoms_added.end());
    std::sort(atoms.begin(), atoms.end());
    graph.atom_nodes.reserve(atoms.size());

    std::unordered_map<AtomBlockId, size_t> atom_local_ids;
    atom_local_ids.reserve(atoms.size());
    for (size_t iatom = 0; iatom < atoms.size(); ++iatom) {
        atom_local_ids.insert({atoms[iatom], iatom + FIRST_ATOM_NODE});
    }

    std::unordered_map<AtomBlockId, int> atom_primitive_nums;
    std::unordered_map<int, int> terminal_to_primitive_num = terminal_primitive_nums(lb_graph);
    for (const t_intra_lb_net& lb_net : intra_lb_nets) {
        VTR_ASSERT(lb_net.atom_pins.size() == lb_net.terminals.size());
        for (size_t ipin = 0; ipin < lb_net.atom_pins.size(); ++ipin) {
            AtomPinId atom_pin = lb_net.atom_pins[ipin];
            if (!atom_pin) {
                continue;
            }

            AtomBlockId atom_id = atom_netlist.pin_block(atom_pin);
            if (atom_local_ids.find(atom_id) == atom_local_ids.end()) {
                continue;
            }

            auto primitive_num_iter = terminal_to_primitive_num.find(lb_net.terminals[ipin]);
            if (primitive_num_iter != terminal_to_primitive_num.end()) {
                atom_primitive_nums.insert({atom_id, primitive_num_iter->second});
            }
        }
    }

    for (AtomBlockId atom_id : atoms) {
        auto primitive_num_iter = atom_primitive_nums.find(atom_id);
        int primitive_num = primitive_num_iter == atom_primitive_nums.end() ? -1 : primitive_num_iter->second;
        graph.atom_nodes.push_back({size_t(atom_id),
                                    models.model_name(atom_netlist.block_model(atom_id)),
                                    primitive_num});
    }

    std::unordered_set<AtomNetId> cluster_net_set;
    for (AtomBlockId atom_id : atoms) {
        for (AtomPinId pin_id : atom_netlist.block_pins(atom_id)) {
            AtomNetId net_id = atom_netlist.pin_net(pin_id);
            if (net_id) {
                cluster_net_set.insert(net_id);
            }
        }
    }

    std::vector<AtomNetId> cluster_nets(cluster_net_set.begin(), cluster_net_set.end());
    std::sort(cluster_nets.begin(), cluster_nets.end());
    for (AtomNetId net_id : cluster_nets) {
        AtomBlockId driver_atom = AtomBlockId::INVALID();
        AtomPinId driver_pin = atom_netlist.net_driver(net_id);
        if (driver_pin) {
            driver_atom = atom_netlist.pin_block(driver_pin);
        }
        auto driver_iter = atom_local_ids.find(driver_atom);
        bool driver_is_internal = driver_iter != atom_local_ids.end();

        for (AtomPinId sink_pin : atom_netlist.net_sinks(net_id)) {
            AtomBlockId sink_atom = atom_netlist.pin_block(sink_pin);
            auto sink_iter = atom_local_ids.find(sink_atom);
            bool sink_is_internal = sink_iter != atom_local_ids.end();

            if (!driver_is_internal && !sink_is_internal) {
                continue;
            }

            size_t source_id = driver_is_internal ? driver_iter->second : EXTERNAL_INPUT_NODE;
            size_t sink_id = sink_is_internal ? sink_iter->second : EXTERNAL_OUTPUT_NODE;
            graph.edges.push_back({source_id, sink_id, size_t(net_id)});
        }
    }

    return graph;
}

/**
 * @brief Streaming writer for the packing-graph JSON array.
 */
class PackingGraphWriter {
  public:
    /** @brief Open a fresh packing_graph.json array. */
    PackingGraphWriter()
        : output_("packing_graph.json", std::ios::out | std::ios::trunc) {
        if (!output_) {
            VTR_LOG_WARN("Failed to open packing_graph.json for writing.\n");
            return;
        }
        output_ << "[\n";
    }

    /** @brief Finalize the JSON array. */
    ~PackingGraphWriter() {
        if (output_) {
            output_ << "\n]\n";
        }
    }

    /** @brief Append one serialized graph record to the output array. */
    void write_record(const t_packing_graph_record& graph) {
        if (!output_) {
            return;
        }
        if (num_records_ != 0) {
            output_ << ",\n";
        }

        output_ << "  {\n"
                << "    \"schema\": \"vtr.intra_lb_route.v1\",\n"
                << "    \"sample_id\": " << num_records_ << ",\n"
                << "    \"cluster_type\": " << json_string(graph.cluster_type) << ",\n"
                << "    \"result\": {\"is_impossible\": " << (graph.is_impossible ? "true" : "false")
                << ", \"route_succeeded\": " << (graph.route_succeeded ? "true" : "false") << "},\n"
                << "    \"nodes\": [\n"
                << "      {\"id\": 0, \"type\": \"external_input\"},\n"
                << "      {\"id\": 1, \"type\": \"external_output\"}";

        for (size_t iatom = 0; iatom < graph.atom_nodes.size(); ++iatom) {
            const t_packing_graph_atom_node& atom = graph.atom_nodes[iatom];
            output_ << ",\n"
                    << "      {\"id\": " << iatom + FIRST_ATOM_NODE
                    << ", \"type\": \"atom\""
                    << ", \"atom_id\": " << atom.atom_id
                    << ", \"model\": " << json_string(atom.model)
                    << ", \"primitive_num\": " << atom.primitive_num << "}";
        }
        output_ << "\n    ],\n"
                << "    \"edges\": [";

        for (size_t iedge = 0; iedge < graph.edges.size(); ++iedge) {
            const t_packing_graph_edge& edge = graph.edges[iedge];
            output_ << (iedge == 0 ? "\n" : ",\n")
                    << "      {\"src\": " << edge.source_node
                    << ", \"dst\": " << edge.sink_node
                    << ", \"atom_net_id\": " << edge.atom_net_id << "}";
        }

        if (!graph.edges.empty()) {
            output_ << '\n';
        }
        output_ << "    ]\n"
                << "  }";
        output_.flush();
        ++num_records_;
    }

  private:
    std::ofstream output_;   ///< Output stream for packing_graph.json.
    size_t num_records_ = 0; ///< Number of records written so far.
};

} // namespace

void write_packing_graph(const std::string& cluster_type,
                         const std::unordered_set<AtomBlockId>& atoms,
                         const std::vector<t_intra_lb_net>& intra_lb_nets,
                         const std::vector<t_lb_type_rr_node>& lb_graph,
                         bool is_impossible,
                         bool route_succeeded,
                         const t_mode_selection_status& mode_status) {
    if (mode_status.is_mode_issue()) {
        return;
    }

    static PackingGraphWriter writer;
    writer.write_record(make_packing_graph(cluster_type,
                                           atoms,
                                           intra_lb_nets,
                                           lb_graph,
                                           is_impossible,
                                           route_succeeded));
}
