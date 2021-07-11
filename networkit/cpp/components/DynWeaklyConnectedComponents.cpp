// no-networkit-format
/*
* DynDynWeaklyConnectedComponents.cpp
*
*  Created on: June 20, 2017
*      Author: Eugenio Angriman
*/

#include <networkit/components/DynWeaklyConnectedComponents.hpp>

namespace NetworKit {

    DynWeaklyConnectedComponents::DynWeaklyConnectedComponents(const Graph& G) : G(&G) {
        if (!G.isDirected()) {
            throw std::runtime_error("Weakly Connected Components can be computeed for directed graphs. Use ConnectedComponents for undirected graphs.");
        }
    }


    void DynWeaklyConnectedComponents::init() {
        edgesMap.clear();
        indexEdges();
        components.assign(G->upperNodeIdBound(), none);
        isTree.assign(edgesMap.size(), false);
        compSize.clear();
        std::queue<index> emptyQueue;
        swap(emptyQueue, componentIds);
        hasRun = false;
    }


    void DynWeaklyConnectedComponents::run() {

        init();

        // Queue for BFS.
        std::queue<node> q;

        // Perform BFSs to assign a component ID to each node.
        G->forNodes([&](node u) {

            // Node u has not been visited.
            if (components[u] == none) {

                // New component ID.
                index c = compSize.size();
                components[u] = c;
                compSize.insert(std::make_pair(c, 1));

                // Start a new BFS from u.
                q.push(u);

                do {
                    node v = q.front();
                    q.pop();

                    // Enqueue neighbors (both from in and out edges) and set
                    // new component.
                    G->forNeighborsOf(v, [&](node w) {
                        updateComponent(c, w, q, v);
                    });

                    G->forInNeighborsOf(v, [&](node w) {
                        updateComponent(c, w, q, v);
                    });
                } while (!q.empty());
            }
        });

        hasRun = true;
    }


    void DynWeaklyConnectedComponents::updateComponent(
        index c,
        node w,
        std::queue<node>& q,
        node v
    ) {
        if (components[w] == none) {
            q.push(w);
            components[w] = c;
            isTree[edgesMap.find(makePair(v, w))->second] = true;
            compSize.find(c)->second += 1;
        }
    }


    void DynWeaklyConnectedComponents::indexEdges() {
        edgeid eid = 0;
        G->forEdges([&] (node u, node v) {
            if (edgesMap.find(makePair(u, v)) == edgesMap.end()) {
                insertEdgeIntoMap(u, v, eid);
                ++eid;
            }
        });
    }


    void DynWeaklyConnectedComponents::insertEdgeIntoMap(
        node u,
        node v,
        edgeid eid
    ) {
        edgesMap.insert(std::make_pair(makePair(u, v), eid));
    }


    edgeid DynWeaklyConnectedComponents::updateMapAfterAddition(
        node u, node v
    ) {

        std::map<std::pair<node, node>, edgeid>::iterator it = edgesMap.find(
            makePair(u, v)
        );


        if (it == edgesMap.end()) {
            // Adding edge never deleted before
            edgeid newId = edgesMap.size();
            insertEdgeIntoMap(u, v, newId);
            return newId;
        }

        return it->second;
    }


    void DynWeaklyConnectedComponents::update(GraphEvent event) {

        assureFinished();

        if (event.type == GraphEvent::EDGE_ADDITION) {
            addEdge(event.u, event.v);
        }
        else if (event.type == GraphEvent::EDGE_REMOVAL) {
            removeEdge(event.u, event.v);
        }
    }


    void DynWeaklyConnectedComponents::updateBatch(
        const std::vector<GraphEvent> &
    ) {
        run();
    }


    void DynWeaklyConnectedComponents::addEdge(node u, node v) {

        edgeid eid = updateMapAfterAddition(u, v);

        // If u and v are already in the same component, we
        // don't have to do anything. Same thing if edge (v, u) already exists.
        index maxComp = std::max(components[u], components[v]);
        index minComp = std::min(components[u], components[v]);

        if (maxComp == minComp || G->hasEdge(v, u)) {
            updateTreeAfterAddition(eid, false);
            return;
        }

        // in the other case, we can merge the two components in an undirected
        // graph
        G->parallelForNodes([&](node w) {
            // We update the component with higher index with the lower index
            if (components[w] == maxComp) {
                components[w] = minComp;
            }
        });

        compSize.find(minComp)->second += compSize.find(maxComp)->second;
        compSize.erase(maxComp);
        componentIds.push(maxComp);

        updateTreeAfterAddition(eid, true);
    }


    void DynWeaklyConnectedComponents::updateTreeAfterAddition(
        edgeid eid, bool partOfTree
    ) {
        if (isTree.size() > eid) {
            isTree[eid] = partOfTree;
        }
        else if (isTree.size() == eid) {
            isTree.push_back(partOfTree);
        }
        else throw std::runtime_error("Edge indexing error");
    }


    void DynWeaklyConnectedComponents::removeEdge(node u, node v) {

        edgeid eid = edgesMap.find(makePair(u, v))->second;

        // If (u, v) is not part of the spanning tree or if edge (v, u) already
        // exists we don't have to do nothing.
        if (!isTree[eid]) {
            return;
        }

        // Edge "eid" is removed from the graph. For performance reasons we
        // keep it in memory, for coherence we claim that it is no more part of
        // the spanning tree.
        isTree[eid] = false;

        index nextId = nextAvailableComponentId(false);

        std::vector<node> newCmp(components);
        newCmp[u] = nextId;
        count newCmpSize = 1;

        tmpDistances.assign(G->upperNodeIdBound(), none);
        std::queue<node> q;
        q.push(u);
        tmpDistances[u] = 0;

        bool connected = false;


        // Berform BFS from v to check if v reaches u
        do {

            node s = q.front();
            q.pop();

            count d = tmpDistances[s] + 1;

            // Enqueue not visited neighbors reachable through outgoing edges
            G->forNeighborsOf(s, [&](node w) {
                if (!connected) {
                    if (tmpDistances[w] == none) {
                        tmpDistances[w] = d;
                        if (w == v) { // Found another path from u to v
                            // Backtracks the path from v to u and marks all its
                            // nodes as part of the spanning tree
                            reverseBFS(u, v);
                            connected = true;
                            return;
                        }

                        newCmp[w] = nextId;
                        ++newCmpSize;
                        q.push(w);
                    }
                }
            });

            // Checking in neighbors (forNeighborsOf gets only the out-degree
            // neighbors).
            if (!connected) {
                G->forInNeighborsOf(s, [&](node w) {
                    if (!connected) {
                        if (tmpDistances[w] == none) {
                            tmpDistances[w] = d;
                            if (w == v) { // Found another path from u to v
                                // Backtracks the path from v to u and marks all
                                // its nodes as part of the spanning tree
                                reverseBFS(u, v);
                                connected = true;
                                return;
                            }

                            newCmp[w] = nextId;
                            ++newCmpSize;
                            q.push(w);
                        }
                    }
                });
            }

            if (connected) {
                break;
            }

        } while (!q.empty());

        if (!connected) {
            // TODO: a more elegant way to assign new ids.
            nextId = nextAvailableComponentId(true);
            compSize.find(components[v])->second -= newCmpSize;
            compSize.insert(std::make_pair(nextId, newCmpSize));
            components = newCmp;
        }
    }

    void DynWeaklyConnectedComponents::reverseBFS(node u, node v) {

        std::queue<node> q;
        q.push(v);
        count d = tmpDistances[v];
        count level = 1;

        do {
            node s = q.front();
            q.pop();

            bool nextEdgeFound = false;
            G->forNeighborsOf(s, [&](node w) {
                if (!nextEdgeFound) {
                    if (visitNodeReversed(
                        u, s, w, v, d, q, nextEdgeFound, level)
                    ) {
                        return;
                    }
                }
            });

            if (!nextEdgeFound) {
                G->forInNeighborsOf(s, [&](node w) {
                    if (!nextEdgeFound) {
                        if (visitNodeReversed(
                            u, s, w, v, d, q, nextEdgeFound, level)
                        ) {
                            return;
                        }
                    }
                });
            }
            ++level;
        } while (!q.empty());
    }


    bool DynWeaklyConnectedComponents::visitNodeReversed(
        node u,
        node s,
        node w,
        node,
        count d,
        std::queue<node>& q,
        bool& nextEdgeFound,
        count level
    ) {

        // Reverse BFS finished
        if (w == u) {
            isTree[edgesMap.find(makePair(w, s))->second] = true;
            nextEdgeFound = true;
            return true;
        }

        // Found next node for reverse BFS
        if ((tmpDistances[w] != none) && (d == tmpDistances[w] + level)) {
            isTree[edgesMap.find(makePair(w, s))->second] = true;
            nextEdgeFound = true;
            q.push(w);
            return true;
        }

        // Discarding node from reverse path
        return false;
    }


    index DynWeaklyConnectedComponents::nextAvailableComponentId(bool eraseId) {
        if (componentIds.empty()) {
            return compSize.size();
        }
        index result = componentIds.front();
        if (eraseId) {
            componentIds.pop();
        }
        return result;
    }


    std::pair<node, node> DynWeaklyConnectedComponents::makePair(
        node u, node v
    ) {
        node from, to;
        if (u > v) {
            from = v;
            to = u;
        }
        else {
            from = u;
            to = v;
        }

        return std::make_pair(from, to);
    }


    std::vector<std::vector<node>> DynWeaklyConnectedComponents::getComponents() {
        assureFinished();

        std::vector<std::vector<node>> result(compSize.size());
        std::map<index, count> compIndex;

        int i = 0;
        for (auto it=compSize.begin(); it!=compSize.end(); ++it) {
            auto indexIterator = compIndex.find(it->first);
            if (indexIterator == compIndex.end()) {
                compIndex.insert(std::make_pair(it->first, i));
                ++i;
            }
        }

        G->forNodes([&](node u) {
            result[compIndex.find(components[u])->second].push_back(u);
        });

        return result;
    }
}
