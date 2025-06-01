#pragma once

#include <unordered_map>
#include <vector>
#include "prepack.h"
// TODO: Add docs

class Prepacker;
class FlatPlacementInfo;
class AtomContext;
class DeviceContext;

struct t_partition_dimension {
  int x;
  int y;
  int z;
};

class NetlistPartition {
  public:
    const std::vector<PackMoleculeId>& molecules(int partition) const;
    int get_partition(PackMoleculeId molecule) const;
    void set_molecule_partition(PackMoleculeId mol, int partition);
    inline t_partition_dimension get_partition_dimensions() const {return partition_dimensions_;}

    NetlistPartition() = delete;
    NetlistPartition(t_partition_dimension partition_dimensions);

  private:
    std::unordered_map<PackMoleculeId, int> partition_map_;
    std::vector<std::vector<PackMoleculeId>> molecules_;
    std::unordered_map<std::string, int> model_count_;
    t_partition_dimension partition_dimensions_;
};

enum class e_partition_type {
  SPATIAL,
  MIN_CUT,
  SPATIAL_MIN_CUT
};

class NetlistPartitioner {
  public:
    NetlistPartitioner(const FlatPlacementInfo& flat_placement_info, const Prepacker& prepacker, const AtomContext& atom_context, const DeviceContext& device_context);
    
    NetlistPartition get_netlist_partition(e_partition_type partition_type, int num_partitions);

  private:
    NetlistPartition get_spatial_partitioning(int num_partitions);
    NetlistPartition get_graph_partitioning(int num_partitions, bool use_placement_info);
    bool should_partition_mol(PackMoleculeId mol_id);

    const Prepacker& prepacker_;
    const FlatPlacementInfo& flat_placement_info_;
    const AtomContext& atom_context_;
    const DeviceContext& device_context_;

    std::unordered_map<std::string, int> model_count_;

};
