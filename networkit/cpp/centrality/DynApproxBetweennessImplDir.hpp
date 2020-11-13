/*
 * DynApproxBetweennssImplDir.hpp
 *
 * Author: Charmaine Ndolo <Charmaine.Ndolo@hu-berlin.de>
 */

// networkit-format
#ifndef CENTRALITY_DYN_APPROX_BETWEENNESS_IMPL_DIR_H
#define CENTRALITY_DYN_APPROX_BETWEENNESS_IMPL_DIR_H

#include <algorithm>
#include <utility>
#include <vector>

#include <networkit/base/DynAlgorithm.hpp>
#include <networkit/centrality/Centrality.hpp>
#include <networkit/distance/DynSSSP.hpp>
#include <networkit/dynamics/GraphEvent.hpp>

namespace NetworKit {

class DynApproxBetweennessImplDir : public Centrality, public DynAlgorithm {

public:
    DynApproxBetweennessImplDir(const Graph &G, double epsilon = 0.1, double delta = 0.1,
                                bool storePredecessors = true);

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

    count computeVDdirected();

private:
    count DFS(node u, count j, Graph &sccDAG, std::vector<bool> &marked, std::vector<node> &sorted);
    void sampleNewPaths(count start, count end);
    bool storePreds = false; // TODO add vector with predecessors
    double epsilon;          //!< maximum error
    double delta;
    const double c = 0.5; // universal positive constant - see reference in paper
    count r;
    count r2;
    std::vector<std::unique_ptr<DynSSSP>> sssp;
    std::vector<std::unique_ptr<DynSSSP>> sssp2;
    std::vector<std::vector<count>> npaths;

    std::vector<count> vis;
    std::vector<node> u;
    std::vector<node> v;
    std::vector<count> maxDist;
    std::vector<count> maxDist2;
    std::vector<std::vector<node>> sampledPaths;

    const count infDist = std::numeric_limits<count>::max();
};
} // namespace NetworKit

#endif
