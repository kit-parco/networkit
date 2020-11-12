/*
 * DynApproxBetweennssImpl.hpp
 *
 * Author: Charmaine Ndolo <Charmaine.Ndolo@hu-berlin.de>
 */

// networkit-format
#ifndef CENTRALITY_DYN_APPROX_BETWEENNESS_IMPL_H
#define CENTRALITY_DYN_APPROX_BETWEENNESS_IMPL_H

#include <algorithm>
#include <utility>
#include <vector>

#include <networkit/base/DynAlgorithm.hpp>
#include <networkit/centrality/Centrality.hpp>
#include <networkit/distance/DynSSSP.hpp>
#include <networkit/dynamics/GraphEvent.hpp>

namespace NetworKit {

class DynApproxBetweennessImpl : public Centrality, public DynAlgorithm {

public:
    DynApproxBetweennessImpl(const Graph &G, double epsilon = 0.1, double delta = 0.1,
                             bool storePredecessors = true, double universalConstant = 1.0);

    /**
     * Runs the static approximated betweenness centrality algorithm on the initial graph.
     */
    void run() override;

    /**
     * Updates the betweenness centralities after an edge insertions on the graph.
     * Notice: it works only with edge insertions and the graph has to be connected.
     *
     * @param e The edge insertions.
     */
    void update(GraphEvent e) override;

    /**
     * Updates the betweenness centralities after a batch of edge insertions on the graph.
     * Notice: it works only with edge insertions and the graph has to be connected.
     *
     * @param batch The batch of edge insertions.
     */
    void updateBatch(const std::vector<GraphEvent> &batch) override;

    /**
     * Get number of path samples used for last calculation
     */
    count getNumberOfSamples();

private:
    bool storePreds = true;
    double epsilon; // maximum error
    double delta;
    double universalConstant;
    count r;
    std::vector<std::unique_ptr<DynSSSP>> sssp;
    std::vector<node> u;
    std::vector<node> v;
    std::vector<std::vector<node>> sampledPaths;
};
} // namespace NetworKit

#endif
