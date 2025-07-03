#pragma once

#include <unordered_map>
#include <vector>
#include "prepack.h"
// TODO: Add docs

class Prepacker;
class FlatPlacementInfo;
class AtomContext;
class DeviceContext;

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
    std::unordered_map<std::string, int> model_count_;
};

enum class e_partition_type {
  NONE,
  SPATIAL
};

class NetlistPartitioner {
  public:
    NetlistPartitioner(const FlatPlacementInfo& flat_placement_info, const Prepacker& prepacker, const AtomContext& atom_context, const DeviceContext& device_context);
    
    NetlistPartition get_netlist_partition(e_partition_type partition_type, int num_partitions);

  private:
    NetlistPartition get_spatial_partitioning(int num_partitions);
    NetlistPartition get_unity_partitioning();
    bool should_partition_mol(PackMoleculeId mol_id);

    const Prepacker& prepacker_;
    const FlatPlacementInfo& flat_placement_info_;
    const AtomContext& atom_context_;
    const DeviceContext& device_context_;

    std::unordered_map<std::string, int> model_count_;
};
