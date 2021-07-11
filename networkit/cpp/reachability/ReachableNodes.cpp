
#include <algorithm>
#include <cassert>
#include <omp.h>

#include <networkit/components/ConnectedComponents.hpp>
#include <networkit/components/StronglyConnectedComponents.hpp>
#include <networkit/reachability/ReachableNodes.hpp>

namespace NetworKit {

ReachableNodes::ReachableNodes(const Graph &G, bool exact) : exact(exact), G(&G) {
    reachableLB.resize(G.upperNodeIdBound());
    if (!exact)
        reachableUB.resize(G.upperNodeIdBound());
}

void ReachableNodes::run() {
    if (G->isDirected())
        runDirected();
    else
        runUndirected();
    hasRun = true;
}

void ReachableNodes::runDirected() {
    StronglyConnectedComponents scc(*G);
    scc.run();

    const count nSCCs = scc.numberOfComponents();
    std::vector<std::vector<bool>> visitedCmpGlobal(omp_get_max_threads());
    visitedCmpGlobal.front().resize(nSCCs);
    std::vector<std::vector<index>> sccVec(nSCCs);

    G->forNodes([&](node u) { sccVec[scc.componentOfNode(u)].push_back(u); });
    Graph sccGraph(nSCCs, false, true);

    // Compute the SCC graph
    for (index cmp = 0; cmp < nSCCs; ++cmp) {
        auto &visitedCmp = visitedCmpGlobal.front();
        for (node u : sccVec[cmp]) {
            G->forNeighborsOf(u, [&](node v) {
                const index neighCmp = scc.componentOfNode(v);

                if (cmp != neighCmp && !visitedCmp[neighCmp]) {
                    visitedCmp[neighCmp] = true;
                    sccGraph.addEdge(cmp, neighCmp);
                }
            });
        }
        sccGraph.forNeighborsOf(cmp, [&](node neighCmp) { visitedCmp[neighCmp] = false; });
    }

    if (exact) {
        for (omp_index i = 1; i < static_cast<omp_index>(omp_get_max_threads()); ++i)
            visitedCmpGlobal[i].resize(nSCCs);

        // BFS from each SCC to count the exact reachable nodes from that SCC
        sccGraph.balancedParallelForNodes([&](node sourceCmp) {
            auto &visitedCmp = visitedCmpGlobal[omp_get_thread_num()];
            std::fill(visitedCmp.begin(), visitedCmp.end(), false);
            visitedCmp[sourceCmp] = true;
            std::queue<node> queue;
            queue.push(sourceCmp);

            count reachableFromSource = 0;
            do {
                const node curCmp = queue.front();
                queue.pop();
                assert(curCmp < sccVec.size());
                reachableFromSource += sccVec[curCmp].size();

                sccGraph.forNeighborsOf(curCmp, [&](node neighCmp) {
                    if (!visitedCmp[neighCmp]) {
                        queue.push(neighCmp);
                        visitedCmp[neighCmp] = true;
                    }
                });

            } while (!queue.empty());

            for (node u : sccVec[sourceCmp])
                reachableLB[u] = reachableFromSource;
        });

    } else {
        // BFS from the biggest SCC, estimate the number of reachable nodes from each SCC.
        const node largestSCC =
            std::max_element(sccVec.begin(), sccVec.end(),
                             [](auto &vec1, auto &vec2) { return vec1.size() < vec2.size(); })
            - sccVec.begin();

        std::vector<count> reachLSCC(nSCCs), reachUSCC(nSCCs);
        auto &reachableFromLargestSCC = visitedCmpGlobal.front();
        std::fill(reachableFromLargestSCC.begin(), reachableFromLargestSCC.end(), false);
        reachableFromLargestSCC[largestSCC] = true;

        std::vector<bool> reachesLargestSCC(nSCCs);
        reachesLargestSCC[largestSCC] = true;

        std::queue<count> queue;
        queue.push(largestSCC);

        do {
            const node curCmp = queue.front();
            queue.pop();

            reachLSCC[largestSCC] += sccVec[curCmp].size();

            sccGraph.forNeighborsOf(curCmp, [&](node neighCmp) {
                if (!reachableFromLargestSCC[neighCmp]) {
                    reachableFromLargestSCC[neighCmp] = true;
                    queue.push(neighCmp);
                }
            });
        } while (!queue.empty());

        reachUSCC[largestSCC] = reachLSCC[largestSCC];
        reachesLargestSCC[largestSCC] = true;

        // So far only the largest SCC has reachUSCC and reachLSCC > 0

        std::vector<count> reachUWithoutLargestSCC(nSCCs);
        // Dynamic programming to compute number of reachable vertices
        sccGraph.forNodes([&](node curCmp) {
            if (curCmp == largestSCC)
                return;

            sccGraph.forNeighborsOf(curCmp, [&](node neighCmp) {
                reachLSCC[curCmp] = std::max(reachLSCC[curCmp], reachLSCC[neighCmp]);

                if (!reachableFromLargestSCC[neighCmp])
                    reachUWithoutLargestSCC[curCmp] += reachUWithoutLargestSCC[neighCmp];

                reachUSCC[curCmp] += reachUSCC[neighCmp];
                reachUSCC[curCmp] = std::min(reachUSCC[curCmp], G->upperNodeIdBound());
                reachesLargestSCC[curCmp] =
                    reachesLargestSCC[curCmp] || reachesLargestSCC[neighCmp];
            });

            if (reachesLargestSCC[curCmp])
                reachUSCC[curCmp] = reachUWithoutLargestSCC[curCmp] + reachUSCC[curCmp];

            reachLSCC[curCmp] += sccVec[curCmp].size();
            reachUSCC[curCmp] += sccVec[curCmp].size();
            reachUSCC[curCmp] = std::min(reachUSCC[curCmp], G->upperNodeIdBound());
        });

        G->parallelForNodes([&](node u) {
            reachableLB[u] = reachLSCC[scc.componentOfNode(u)];
            reachableUB[u] = reachUSCC[scc.componentOfNode(u)];
        });
    }
}

void ReachableNodes::runUndirected() {
    ConnectedComponents cc(*G);
    cc.run();
    const auto componentSizes = cc.getComponentSizes();

    G->parallelForNodes([&](node u) { reachableLB[u] = componentSizes.at(cc.componentOfNode(u)); });
}
} // namespace NetworKit
