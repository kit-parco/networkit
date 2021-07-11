/*
 * IncompleteDijkstra.h
 *
 *  Created on: 15.07.2014
 *      Author: dhoske
 */

#ifndef NETWORKIT_DISTANCE_INCOMPLETE_DIJKSTRA_HPP_
#define NETWORKIT_DISTANCE_INCOMPLETE_DIJKSTRA_HPP_

#include <unordered_set>
#include <vector>

#include <tlx/container/d_ary_addressable_int_heap.hpp>

#include <networkit/auxiliary/PrioQueue.hpp>
#include <networkit/distance/IncompleteSSSP.hpp>
#include <networkit/graph/Graph.hpp>

namespace NetworKit {

/**
 * @ingroup distance
 * Implementation of @a IncompleteSSSP using a normal
 * Dijkstra with binary heaps.
 */
class IncompleteDijkstra : public IncompleteSSSP {
public:
    /**
     * Creates a IncompleteDijkstra instance from the sources in
     * @a sources (act like a super source) in the graph @a G.
     * The edges in @a G must have nonnegative weight and @a G should
     * not be null.
     *
     * We also consider the nodes in @a explored to not exist
     * if @a explored is not null.
     *
     * @warning We do not copy @a G or @a explored, but store a
     * non-owning pointer to them. Otherwise IncompleteDijkstra would not
     * be more efficient than normal Dijkstra. Thus, @a G and @a explored
     * must exist at least as long as this IncompleteDijkstra instance.
     *
     * @todo This is somewhat ugly, but we do not want introduce a
     * std::shared_ptr<> since @a G and @a explored could well
     * be stack allocated.
     */
    IncompleteDijkstra(const Graph *G, const std::vector<node> &sources,
                       const std::unordered_set<node> *explored = nullptr);

    bool hasNext() override;
    std::pair<node, edgeweight> next() override;

private:
    // Stored reference to outside data structures
    const Graph *G;
    const std::unordered_set<node> *explored;

    std::vector<edgeweight> dists;

    struct CompareDistance {
        CompareDistance(const std::vector<edgeweight> *distance) : distance(distance) {}
        bool operator()(node x, node y) const noexcept { return (*distance)[x] < (*distance)[y]; }

    private:
        const std::vector<edgeweight> *distance;
    };

    tlx::d_ary_addressable_int_heap<node, 2, CompareDistance> heap;
};

} // namespace NetworKit

#endif // NETWORKIT_DISTANCE_INCOMPLETE_DIJKSTRA_HPP_
