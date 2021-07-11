#ifndef NETWORKIT_SCD_SET_CONDUCTANCE_HPP_
#define NETWORKIT_SCD_SET_CONDUCTANCE_HPP_

#include <set>

#include <networkit/base/Algorithm.hpp>
#include <networkit/graph/Graph.hpp>

namespace NetworKit {

/**
 * Calculates the conductance of a set of nodes, i.e., the weight of all edges
 * between the set and the rest of the graph divided by the minimum
 * of the volume (the sum of the weighted degrees) of the community and the rest
 * of the graph.
 */
class SetConductance : public Algorithm {
public:
    /**
     * Construct the SetConductance with the given graph and community.
     *
     * @param G The graph
     * @param community The set of nodes to examine.
     */
    SetConductance(const Graph &g, const std::set<node> &community);

    /**
     * Calculate the conductance.
     */
    void run() override;

    /**
     * Get the calculated conductance score.
     *
     * @return The conductance.
     */
    double getConductance() const;

private:
    const Graph *g;
    const std::set<node> *community;
    double conductance;
};

} // namespace NetworKit

#endif // NETWORKIT_SCD_SET_CONDUCTANCE_HPP_
