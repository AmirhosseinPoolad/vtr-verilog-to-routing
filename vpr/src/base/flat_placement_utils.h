/**
 * @file
 * @author  Alex Singer
 * @date    March 2025
 * @brief   Utility methods for working with flat placements.
 */

#pragma once

#include <cstdlib>
#include "flat_placement_types.h"

/**
 * @brief Returns the manhattan distance (L1 distance) between two flat
 *        placement locations.
 */
inline float get_manhattan_distance(const t_flat_pl_loc& loc_a,
                                    const t_flat_pl_loc& loc_b) {
    return std::abs(loc_a.x - loc_b.x) + std::abs(loc_a.y - loc_b.y) + std::abs(loc_a.layer - loc_b.layer);
}

/**
 * @brief Returns the midpoint between two flat
 *        placement locations.
 */
inline t_flat_pl_loc get_midpoint(const t_flat_pl_loc& loc_a,
                           const t_flat_pl_loc& loc_b) {
    return t_flat_pl_loc({(loc_a.x + loc_b.x)/2, (loc_a.y + loc_b.y)/2, (loc_a.layer + loc_b.layer)/2});
}
