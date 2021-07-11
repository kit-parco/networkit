// no-networkit-format
/*
 * UnionFind.hpp
 *
 *  Created on: 11.08.2014
 *      Author: Marcel Radermacher
 *      Changed a bit by Henning Meyerhenke to reflect union by rank and path compression
 *        as taught in "Algorithms 1"
 */

#ifndef NETWORKIT_STRUCTURES_UNION_FIND_HPP_
#define NETWORKIT_STRUCTURES_UNION_FIND_HPP_

#include <vector>

#include <networkit/Globals.hpp>
#include <networkit/structures/Partition.hpp>

namespace NetworKit {

/**
 * @ingroup structures
 * Implements the Union Find data structure to maintain disjoint sets efficiently.
 * Uses path compression and union by rank to achieve running time linear in
 * the number of elements times the inverse Ackermann function.
 */
class UnionFind final {
    std::vector<index> parent;
    std::vector<unsigned char> rank;

public:

    /**
     * Create a new set representation with not more than @p max_element elements.
     * Initially every element is in its own set.
     * @param max_element maximum number of elements
     */
    UnionFind(index max_element) : parent(max_element), rank(max_element, 0) {
        allToSingletons();
    }

    /**
     * Assigns every element to a singleton set.
     * Set id is equal to element id.
     */

    void allToSingletons();

    /**
     * Find the representative to element @u
     * @param u element
     * @return representative of set containing @u
     */
    index find(index u);

    /**
     *  Merge the two sets contain @u and @v
     *  @param u element u
     *  @param v element v
     */
    void merge(index u, index v);

    /**
     * Convert the Union Find data structure to a Partition
     * @return Partition equivalent to the union find data structure
     * */
    Partition toPartition();
};
} /* namespace NetworKit */
#endif // NETWORKIT_STRUCTURES_UNION_FIND_HPP_
