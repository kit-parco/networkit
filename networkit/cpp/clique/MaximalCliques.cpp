// no-networkit-format
#include <algorithm>
#include <cassert>

#include <utility>
#include <networkit/auxiliary/SignalHandling.hpp>
#include <networkit/centrality/CoreDecomposition.hpp>
#include <networkit/clique/MaximalCliques.hpp>

namespace {
    // Private implementation namespace
    using NetworKit::node;
    using NetworKit::count;
    using NetworKit::index;

    class MaximalCliquesImpl {
    private:
        const NetworKit::Graph* G;
        std::vector<std::vector<node>>* result;
        std::function<void(const std::vector<node>&)>& callback;
        bool maximumOnly;
        count maxFound;

        std::vector<node> pxvector;
        std::vector<node> pxlookup;

        std::vector<index> firstOut;
        std::vector<node> head;

    public:
        MaximalCliquesImpl(const NetworKit::Graph& G, std::vector<std::vector<node>>& result,
                std::function<void(const std::vector<node>&)>& callback, bool maximumOnly) :
            G(&G), result(&result), callback(callback), maximumOnly(maximumOnly), maxFound(0),
            pxvector(G.numberOfNodes()), pxlookup(G.upperNodeIdBound()),
            firstOut(G.upperNodeIdBound() + 1), head(G.numberOfEdges()) {}

    private:
        void buildOutGraph() {
            index currentOut = 0;
            for (node u = 0; u < G->upperNodeIdBound(); ++u) {
                firstOut[u] = currentOut;
                if (G->hasNode(u)) {
                    index xpboundU = pxlookup[u];
                    G->forEdgesOf(u, [&](node v) {
                        if (xpboundU < pxlookup[v]) {
                            head[currentOut++] = v;
                        }
                    });
                }
            }
            firstOut[G->upperNodeIdBound()] = currentOut;
        }

        template <typename F>
        void forOutEdgesOf(node u, F callback) const {
            for (index i = firstOut[u]; i < firstOut[u + 1]; ++i) {
                callback(head[i]);
            }
        }

        bool hasNeighbor(node u, node v) const {
            for (index i = firstOut[u]; i < firstOut[u + 1]; ++i) {
                if (head[i] == v) return true;
            }

            return false;
        }

        count outDegree(node u) const {
            return firstOut[u + 1] - firstOut[u];
        }

        void swapNodeToPos(node u, index pos) {
            assert(pos < pxvector.size());
            node pxvec2 = pxvector[pos];
            std::swap(pxvector[pxlookup[u]], pxvector[pos]);
            pxlookup[pxvec2] = pxlookup[u];
            pxlookup[u] = pos;
        }

    public:
        void run() {
            NetworKit::CoreDecomposition cores(*G, false, false, true);
            cores.run();

            Aux::SignalHandler handler;
            handler.assureRunning();

            const auto& orderedNodes = cores.getNodeOrder();

            index ii = 0;
            for (const node u : orderedNodes) {
                pxvector[ii] = u;
                pxlookup[u] = ii;
                ii += 1;
            }

            handler.assureRunning();

#ifndef NDEBUG
            for (auto u : orderedNodes) {
                assert(pxvector[pxlookup[u]] == u);
            }
#endif

            // Store out-going neighbors in the direction of higher core numbers.
            // This means that the out-degree is bounded by the maximum core number.
            buildOutGraph();

            handler.assureRunning();

            for (index iu = orderedNodes.size(); iu-- > 0; ) {
                node u = orderedNodes[iu];
                index xpbound = iu + 1;

#ifndef NDEBUG
                for (auto v : orderedNodes) {
                    if (v == u)
                        break;

                    assert(pxlookup[v] < xpbound);
                }
#endif

                // Check if u can be the starting point of a new clique
                // of size greater than maxFound.
                // Note that the clique starting at u could be of
                // size outDegree(u) + 1, but then it is still only the
                // same size as maxFound.
                if (maximumOnly && maxFound > outDegree(u)) {
                    swapNodeToPos(u, iu);
                    continue;
                }

                count xcount = 0;
                count pcount = 0;
                G->forNeighborsOf(u, [&] (node v) {

#ifndef NDEBUG
                    assert(pxlookup[v] < pxvector.size());

                    assert(xcount <= xpbound);
                    assert(pcount <= pxvector.size() - xpbound);
#endif

                    if (pxlookup[v] < xpbound) { // v is in X
                        swapNodeToPos(v, xpbound - xcount - 1);
                        xcount += 1;
                    } else { // v is in P
                        swapNodeToPos(v, xpbound + pcount);
                        pcount += 1;
                    }
                });

#ifndef NDEBUG
                { // assert all neighbors of u were stored in one range around xpbound

                    assert(xcount + pcount == G->degree(u));

                    std::vector<index> neighborPositions;
                    G->forNeighborsOf(u, [&](node v) {
                        neighborPositions.push_back(pxlookup[v]);
                        assert(pxvector[pxlookup[v]] == v);
                    });

                    std::sort(neighborPositions.begin(), neighborPositions.end());
                    for (index i = 0; i < neighborPositions.size(); ++i) {
                        assert(neighborPositions[i] == xpbound - xcount + i);
                    }
                }
#endif

                std::vector<node> r = {u};
                tomita(xpbound - xcount, xpbound, xpbound + pcount, r);

                swapNodeToPos(u, iu);

                handler.assureRunning();
            }
        }

        void tomita(index xbound, index xpbound, index pbound, std::vector<node>& r) {
            if (xbound == pbound) { //if (X, P are empty)
                if (callback) {
                    callback(r);
                } else if (!maximumOnly) {
                    result->push_back(r);
                } else if (r.size() > maxFound) {
                    result->clear();
                    result->push_back(r);
                    maxFound = r.size();
                }
                return;
            }

            if (xpbound == pbound) return;

#ifndef NDEBUG
            assert(xbound <= xpbound);
            assert(xpbound <= pbound);
            assert(pbound <= pxvector.size());
#endif

            Aux::SignalHandler handler;
            handler.assureRunning();

            node u = findPivot(xbound, xpbound, pbound);
            std::vector<node> movedNodes;

            // Find all nodes in P that are not neighbors of the pivot
            // this step is necessary as the next loop changes pxvector,
            // which prohibits iterating over it in the same loop.
            std::vector<node> toCheck;

            // Step 1: mark all outgoing neighbors of the pivot in P
            std::vector<bool> pivotNeighbors(pbound - xpbound);
            forOutEdgesOf(u, [&](node v) {
                index vpos = pxlookup[v];
                if (vpos >= xpbound && vpos < pbound) {
                    pivotNeighbors[vpos - xpbound] = true;
                }
            });

            // Step 2: for all not-yet marked notes check if they have the pivot as neighbor.
            // If not: they are definitely a non-neighbor.
            for (index i = xpbound; i < pbound; i++) {
                if (!pivotNeighbors[i - xpbound]) {
                    node p = pxvector[i];

                    if (!hasNeighbor(p, u)) {
                        toCheck.push_back(p);
                    }
                }
            }

            for (auto pxveci : toCheck) {
                count xcount = 0, pcount = 0;

                // Group all neighbors of pxveci in P \cup X around xpbound.
                // Step 1: collect all outgoing neighbors of pxveci
                forOutEdgesOf(pxveci, [&](node v) {
                    if (pxlookup[v] < xpbound && pxlookup[v] >= xbound) { // v is in X
                        swapNodeToPos(v, xpbound - xcount - 1);
                        xcount += 1;
                    } else if (pxlookup[v] >= xpbound && pxlookup[v] < pbound){ // v is in P
                        swapNodeToPos(v, xpbound + pcount);
                        pcount += 1;
                    }
                });

                // Step 2: collect all nodes in X that have not yet been collected
                // and that have pxveci as outgoing neighbor.
                for (index i = xbound; i < xpbound;) {
                    // stop if we have reached the collected neighbors
                    if (i == xpbound - xcount) break;
                    node x = pxvector[i];

                    if (hasNeighbor(x, pxveci)) {
                        swapNodeToPos(x, xpbound - xcount - 1);
                        xcount += 1;
                    } else {
                        // Advance only if we did not swap otherwise we have already
                        // a next candidate at position i.
                        ++i;
                    }
                }

                // Step 3: collect all nodes in P that have not yet been collected
                // and that have pxveci as outgoing neighbor.
                for (index i = xpbound + pcount; i < pbound; ++i) {
                    node p = pxvector[i];

                    if (hasNeighbor(p, pxveci)) {
                        swapNodeToPos(p, xpbound + pcount);
                        pcount += 1;
                    }
                }

                r.push_back(pxveci);

#ifndef NDEBUG
                assert(xpbound + pcount <= pbound);
                assert(xpbound - xcount >= xbound);
#endif

                // only the pcount nodes in P are candidates for the clique,
                // therefore r.size() + pcount is an upper bound for the maximum
                // size of the clique that can still be found in this branch
                // of the recursion.
                if (!maximumOnly || maxFound < (r.size() + pcount)) {
                    tomita(xpbound - xcount, xpbound, xpbound + pcount, r);
                }

                r.pop_back();

                swapNodeToPos(pxveci, xpbound);
                xpbound += 1;
                assert(pxvector[xpbound - 1] == pxveci);
                movedNodes.push_back(pxveci);
            }

            for (node v : movedNodes) {
                //move from X -> P
                swapNodeToPos(v, xpbound - 1);
                xpbound -= 1;
            }

#ifndef NDEBUG
            for (node v : movedNodes) {
                assert(pxlookup[v] >= xpbound);
                assert(pxlookup[v] < pbound);
            }
#endif
        }

        node findPivot(index xbound, index xpbound, index pbound) const {
            // Counts for every node in X \cup P how many outgoing neighbors it has in P
            std::vector<count> pivotNeighbors(pbound - xbound);
            const count psize = pbound-xpbound;

            // Step 1: for all nodes in X count how many outgoing neighbors they have in P
            for (index i = 0; i < xpbound - xbound; i++) {
                node u = pxvector[i + xbound];
                forOutEdgesOf(u, [&](node v) {
                    if (pxlookup[v] >= xpbound && pxlookup[v] < pbound) {
                        ++pivotNeighbors[i];
                    }
                });

                // If a node has |P| neighbors, we cannot find a better candidate
                if (pivotNeighbors[i] == psize) return u;
            }

            // Step 2: for all nodes in P
            // a) increase counts for every neighbor in P \cup X to account for incoming neighbors
            // b) count all outgoing neighbors in P
            for (index i = xpbound - xbound; i < pivotNeighbors.size(); ++i) {
                node u = pxvector[i + xbound];
                forOutEdgesOf(u, [&](node v) {
                    index neighborPos = pxlookup[v];
                    if (neighborPos >= xbound && neighborPos < pbound) {
                        ++pivotNeighbors[neighborPos-xbound];

                        if (neighborPos >= xpbound) {
                            ++pivotNeighbors[i];
                        }
                    }
                });
            }

            node maxnode = pxvector[xbound];
            count maxval = pivotNeighbors[0];

            // Step 3: find maximum
            for (index i = 1; i < pivotNeighbors.size(); ++i) {
                if (pivotNeighbors[i] > maxval) {
                    maxval = pivotNeighbors[i];
                    maxnode = pxvector[i + xbound];
                }
            }

#ifndef NDEBUG
            assert(maxnode < G->upperNodeIdBound() + 1);
#endif

            return maxnode;
        }

    };

}

namespace NetworKit {

MaximalCliques::MaximalCliques(const Graph& G, bool maximumOnly) : G(&G), maximumOnly(maximumOnly) {
}

MaximalCliques::MaximalCliques(const Graph &G,
                               std::function<void(const std::vector<node> &)> callback)
    : G(&G), callback(std::move(callback)), maximumOnly(false) {}

const std::vector<std::vector<node>>& MaximalCliques::getCliques() const {
    if (callback) throw std::runtime_error("MaximalCliques used with callback does not store cliques");
    assureFinished();
    return result;
}

void MaximalCliques::run() {
    hasRun = false;

    result.clear();

    MaximalCliquesImpl(*G, result, callback, maximumOnly).run();

    hasRun = true;
}

}
