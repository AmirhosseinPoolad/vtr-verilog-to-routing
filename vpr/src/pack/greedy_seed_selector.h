#pragma once
/**
 * @file
 * @author  Alex Singer
 * @date    November 2024
 * @brief   The declaration of the greedy seed selector class which selects the
 *          seed molecules for starting new clusters in the greedy clusterer.
 */

#include "prepack.h"
#include "vpr_types.h"
#include "vtr_vector.h"

// Forward declarations
class AtomNetlist;
class ClusterLegalizer;
class LogicalModels;
class PreClusterTimingManager;
class RamMapper;
struct t_molecule_stats;

/**
 * @brief Precomputes eigenvector centrality for the atoms in a netlist.
 *
 * Each driver-to-sink pin connection adds an undirected edge of unit weight.
 * Parallel connections accumulate; self-connections are excluded. All nets,
 * including clock nets, participate. Scores are the non-negative, L2-normalized
 * dominant eigenvector of this adjacency matrix. Isolated atoms score zero.
 * Disconnected components compete for the global dominant eigenvalue; equal
 * dominant eigenvalues are resolved by a uniform initial vector.
 */
class AtomEigenvectorCentrality {
  public:
    /** @brief Compute scores using shifted sparse power iteration with Eigen.
     * Uses a relative residual tolerance of 1e-6 and at most 1000 iterations.
     * Warns and retains the approximation if the iteration limit is reached.
     */
    explicit AtomEigenvectorCentrality(const AtomNetlist& atom_netlist);

    /// @brief Return the precomputed score for a valid atom in the input netlist.
    float get_centrality(AtomBlockId atom_id) const;

  private:
    vtr::vector<AtomBlockId, float> centrality_; ///< Atom scores; invalid ID slots contain -1.
};

/**
 * @brief A selector class which will propose good seed values to start new
 *        clusters in the greedy clusterer.
 *
 * In greedy clustering algorithms, the order in which clusters are generated
 * can have an effect on the quality of the clustering. This class proposes
 * good seed molecules based on heuristics which give each molecule a "seed
 * gain". This class will not propose a molecule which it has already proposed
 * or has already been clustered.
 */
class GreedySeedSelector {
  public:
    /**
     * @brief Constructor of the Greedy Seed Selector class. Pre-computes the
     *        gains of each molecule internally to make getting seeds later very
     *        quick.
     *
     *  @param atom_netlist
     *              The netlist of atoms to cluster.
     *  @param prepacker
     *              The prepacker used to generate pack-pattern molecules of the
     *              atoms in the netlist.
     *  @param seed_type
     *              Controls the algorithm used to compute the seed gain for
     *              each molecule.
     *  @param max_molecule_stats
     *              The maximum stats over all molecules. Used for normalizing
     *              terms in the gain. The local copy is augmented with the
     *              maximum atom eigenvector centrality for BLEND2.
     *  @param pre_cluster_timing_manager
     *              Timing manager class for the primitive netlist. Used to
     *              compute the criticalities of seeds.
     *  @param ram_mapper
     *              The RAM mapper which contains the pre-assigned RAM groups.
     *              If there are any RAM groups, RAM seeds are moved to the
     *              front of the seed list so that RAM clusters are formed first.
     */
    GreedySeedSelector(const AtomNetlist& atom_netlist,
                       const Prepacker& prepacker,
                       const e_cluster_seed seed_type,
                       t_molecule_stats max_molecule_stats,
                       const LogicalModels& models,
                       const PreClusterTimingManager& pre_cluster_timing_manager,
                       const RamMapper& ram_mapper);

    /**
     * @brief Propose a new seed molecule to start a new cluster with. If no
     *        unclustered molecules exist, will return an invalid ID.
     *
     * This method will never return a molecule which has already been clustered
     * (according to the cluster legalizer) and will never propose a molecule
     * that it already proposed.
     *
     * This method assumes that once a molecule is clustered, it will never be
     * unclustered.
     *
     *  @param cluster_legalizer
     *              The cluster legalizer object that is used to create the
     *              clusters. This is used to check if a molecule has already
     *              been clustered or not.
     */
    PackMoleculeId get_next_seed(const ClusterLegalizer& cluster_legalizer);

    // TODO: Maybe create an update_seed_gains method to update the seed molecules
    //       list using current clustering information.

  private:
    /// @brief The index of the next seed to propose in the seed_mols_ vector.
    ///        This is set to 0 in the constructor and incremented as more seeds
    ///        are proposed.
    size_t seed_index_;

    /// @brief A list of seed molecules, sorted in decreasing order of gain. This
    ///        is computed in the constructor and is traversed when a new seed
    ///        is being proposed.
    std::vector<PackMoleculeId> seed_mols_;
};
