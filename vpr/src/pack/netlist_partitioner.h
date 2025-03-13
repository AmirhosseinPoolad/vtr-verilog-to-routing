#pragma once

#include <unordered_map>
#include <vector>
#include "prepack.h"
// TODO: Add docs

class Prepacker;
class FlatPlacementInfo;
class AtomContext;

class NetlistPartition {
  public:
    const std::vector<PackMoleculeId>& molecules(int partition) const;
    int get_partition(PackMoleculeId molecule) const;
    void set_molecule_partition(PackMoleculeId mol, int partition);
    NetlistPartition() = delete;
    NetlistPartition(int num_partitions);

  private:
    std::unordered_map<PackMoleculeId, int> partition_map_;
    std::vector<std::vector<PackMoleculeId>> molecules_;
};

// typedef std::unordered_map<PackMoleculeId, int> NetlistPartition;

enum class e_partition_type {
  SPATIAL,
  MIN_CUT,
  SPATIAL_MIN_CUT
};

class NetlistPartitioner {
  public:
    NetlistPartitioner(const FlatPlacementInfo& flat_placement_info, const Prepacker& prepacker, const AtomContext& atom_context);
    
    NetlistPartition get_netlist_partition(e_partition_type partition_type, int num_partitions);

  private:
    NetlistPartition get_spatial_partitioning(int num_partitions);
    NetlistPartition get_graph_partitioning(int num_partitions, bool use_placement_info);
    const Prepacker& prepacker_;
    const FlatPlacementInfo& flat_placement_info_;
    const AtomContext& atom_context_;
};
