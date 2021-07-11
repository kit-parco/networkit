#ifndef NETWORKIT_SCD_LOCALDEGREEDIRECTEDGRAPH_HPP
#define NETWORKIT_SCD_LOCALDEGREEDIRECTEDGRAPH_HPP

#include <unordered_map>
#include <vector>

#include <networkit/Globals.hpp>
#include <networkit/graph/Graph.hpp>

namespace NetworKit {

/**
 * Graph for local community detection that stores edges in directed
 * form in just one direction depending on the node's degrees and
 * thus allows efficient triangle listing.
 */
template <bool IsWeighted, typename NodeAddedCallbackType>
class LocalDegreeDirectedGraph {
protected:
    // local to global (input graph) id mapping
    std::vector<node> localToGlobalId;
    // map from input graph to local id
    std::unordered_map<node, node> globalToLocalId;

    // for the local graph: outgoing neighbors
    std::vector<node> head;
    // for the local graph: weight of outgoing edges, only used for weighted graphs (otherwise
    // empty)
    std::vector<double> headWeight;

    struct LocalNode {
        /** first index in head where neighbors of the node are stored */
        size_t firstHead;
        /** index after the last index in head where neighbors of the node are stored */
        size_t lastHead;
    };

    // stores for every local node the information where the outgoing neighbors begins and ends
    std::vector<LocalNode> headInfo;
    // stores the degree of every local node in the local graph
    std::vector<count> degree;

    std::vector<node> currentLocalNeighborIds;

    const Graph &g;

    // data structures only for neighbors of a node
    std::vector<bool> isNeighbor;
    std::vector<double> neighborWeight;

    // callback to call whenever a node is added
    NodeAddedCallbackType nodeAddedCallback;

public:
    /**
     * Initialize the local graph as empty graph.
     *
     * @param g The graph from which nodes can be added
     * @param node_added_callback Callback that is called whenever a node is added
     */
    LocalDegreeDirectedGraph(const Graph &g, NodeAddedCallbackType nodeAddedCallback)
        : g(g), nodeAddedCallback(nodeAddedCallback) {
        if (g.isDirected()) {
            throw std::runtime_error("Directed graphs are not supported");
        }
    }

    ~LocalDegreeDirectedGraph() = default;

    /**
     * Map a local node id to a global one
     */
    node toGlobal(node lu) const { return localToGlobalId[lu]; }

    /**
     * Check if the given node exists.
     */
    bool hasNode(node u) const { return globalToLocalId.find(u) != globalToLocalId.end(); }

    /**
     * Ensure that the given global node id exists in the local graph
     *
     * @param u The global node id to check
     * @return The local node id
     */
    node ensureNodeExists(node u) {
        auto gtlIt = globalToLocalId.find(u);

        if (gtlIt == globalToLocalId.end()) {
            node localId = localToGlobalId.size();
            localToGlobalId.push_back(u);
            globalToLocalId[u] = localId;

            double weightedDegree = 0;

            size_t myDegree = 0;

            index myBegin = head.size();
            index myEnd = myBegin;
            index nh = myEnd; // next head if all potential out neighbors were inserted

            g.forEdgesOf(u, [&](node, node v, edgeweight weight) {
                auto it = globalToLocalId.find(v);

                if (it != globalToLocalId.end()) {
                    ++myDegree;

                    auto &nDegree = degree[it->second];
                    ++nDegree;

                    head.push_back(it->second);
                    if (IsWeighted) {
                        headWeight.push_back(weight);
                    }

                    ++myEnd;
                }

                weightedDegree += weight;
                ++nh;
            });

            assert(myEnd == head.size());
            assert(!IsWeighted || myEnd == headWeight.size());
            assert(nh >= myEnd);

            head.resize(nh);
            if (IsWeighted) {
                headWeight.resize(nh);
            }

            degree.push_back(myDegree);

            // Iterate over all neighbors in the local graph
            for (index i = myBegin; i < myEnd;) {
                node ln = head[i];
                edgeweight weight = IsWeighted ? headWeight[i] : defaultEdgeWeight;

                auto &nDegree = degree[ln];
                auto &nInfo = headInfo[ln];

                // Before into the neighbors of ln, check if it can get rid of any
                // neighbors due to being a higher-degree node now.
                for (index ni = nInfo.firstHead; ni < nInfo.lastHead;) {
                    node nn = head[ni];
                    if (nDegree <= degree[nn]) {
                        ++ni;
                    } else {
                        // Move edge to node nn
                        // Store last neighbor of ln at position ni instead, i.e., remove nn
                        head[ni] = head[--nInfo.lastHead];

                        if (IsWeighted) {
                            // Copy weight to neighbor list of nn
                            headWeight[headInfo[nn].lastHead] = headWeight[ni];
                            // Copy weight of last neighbor of nn at position ni
                            headWeight[ni] = headWeight[nInfo.lastHead];
                        }

                        // Store ln as neighbor at the last position of the neighbor list of nn
                        head[headInfo[nn].lastHead++] = ln;
                    }
                }

                if (myDegree < nDegree) {
                    ++i;
                } else {
                    // Add local_id to ln as neighbor instead of adding ln as neighbor here
                    if (IsWeighted) {
                        headWeight[nInfo.lastHead] = weight;
                    }

                    head[nInfo.lastHead++] = localId;

                    // This removes one neighbor
                    --myEnd;
                    // Move the last neighbor at the current position
                    head[i] = head[myEnd];
                    if (IsWeighted) {
                        headWeight[i] = headWeight[myEnd];
                    }

                    assert(ln + 1 == localId
                           || headInfo[ln].lastHead <= headInfo[ln + 1].firstHead);
                }
            }

            headInfo.emplace_back(LocalNode{myBegin, myEnd});

            isNeighbor.push_back(false);

            if (IsWeighted)
                neighborWeight.push_back(0);

            nodeAddedCallback(u, localId, weightedDegree);

            validateEdgesExist();

            return localId;
        } else {
            return gtlIt->second;
        }
    }

    /**
     * Iterate over all triangles that include the global node @a u
     *
     * @param u The node from which triangles shall be listed
     * @param callback The callback to call. Must take two node ids of the triangle as well as three
     * edge weights as parameters.
     */
    template <typename F>
    void forTrianglesOf(node u, F callback) {
        forTrianglesOf(
            u, [](node, edgeweight) {}, callback);
    }

    /**
     * Iterate over all triangles that include the global node @a u
     *
     * @param u The node from which triangles shall be listed
     * @param initNeighborsCallback The callback to call for every local neighbor (local id) as they
     * are collected
     * @param triangleCallback The callback to call for each
     * triangle. Must take two node ids of the triangle as well as
     * three edge weights as parameters.
     */
    template <typename Fi, typename Ft>
    void forTrianglesOf(node u, Fi initNeighborsCallback, Ft triangleCallback) {
        if (g.degree(u) == 0)
            return;

        currentLocalNeighborIds.clear();
        currentLocalNeighborIds.reserve(g.degree(u));

        g.forNeighborsOf(u, [&](node, node v, edgeweight weight) {
            node localId = ensureNodeExists(v);
            currentLocalNeighborIds.emplace_back(localId);
            isNeighbor[localId] = true;

            if (IsWeighted) {
                neighborWeight[localId] = weight;
            }

            initNeighborsCallback(localId, weight);
        });

        for (node lv : currentLocalNeighborIds) {
            for (index ti = headInfo[lv].firstHead; ti < headInfo[lv].lastHead; ++ti) {
                node y = head[ti];

                if (isNeighbor[y]) {
                    edgeweight weightUv = defaultEdgeWeight;
                    edgeweight weightUy = defaultEdgeWeight;
                    edgeweight weightVy = defaultEdgeWeight;
                    if (IsWeighted) {
                        weightUv = neighborWeight[lv];
                        weightUy = neighborWeight[y];
                        weightVy = headWeight[ti];
                    }

                    triangleCallback(lv, y, weightUv, weightUy, weightVy);
                }
            }
        }

        for (node v : currentLocalNeighborIds) {
            isNeighbor[v] = false;
        }
    }

    void validateEdgesExist() {
#ifdef NETWORKIT_SANITY_CHECKS
#ifndef NDEBUG
        g.forNodes([&](node gu) {
            auto gtlIt = globalToLocalId.find(gu);
            if (gtlIt != globalToLocalId.end()) {
                node lu = gtlIt->second;

                g.forNeighborsOf(gu, [&](node gv) {
                    auto gtlV = globalToLocalId.find(gv);
                    if (gtlV != globalToLocalId.end()) {
                        node lv = gtlV->second;

                        bool foundU = false;
                        for (size_t i = headInfo[lu].firstHead; i < headInfo[lu].lastHead; ++i) {
                            if (head[i] == lv) {
                                foundU = true;
                                break;
                            }
                        }
                        bool foundV = false;
                        for (size_t i = headInfo[lv].firstHead; i < headInfo[lv].lastHead; ++i) {
                            if (head[i] == lu) {
                                foundV = true;
                                break;
                            }
                        }

                        assert(foundU != foundV);
                    }
                });
            }
        });
#endif
#endif
    }

    /**
     * Iterate over local neighbors of the node previously passed to forTrianglesOf
     *
     * @param callback The call back that shall be called for each neighbor, needs to take a node
     * and an edgeweight.
     */
    template <typename F>
    void forLocalNeighbors(F callback) {
        for (node v : currentLocalNeighborIds) {
            if (IsWeighted) {
                callback(v, neighborWeight[v]);
            } else {
                callback(v, defaultEdgeWeight);
            }
        }
    }
};

} // namespace NetworKit
#endif
