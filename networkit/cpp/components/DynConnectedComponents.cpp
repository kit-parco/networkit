// no-networkit-format
/*
* DynConnectedComponents.cpp
*
*  Created on: June 2017
*      Author: Eugenio Angriman
*/

#include <networkit/components/DynConnectedComponents.hpp>

namespace NetworKit {

    DynConnectedComponents::DynConnectedComponents(const Graph& G) :
    G(&G) {
        if (G.isDirected()) {
            throw std::runtime_error("Error, connected components of directed graphs cannot be computed, use StronglyConnectedComponents instead.");
        }
    }


    void DynConnectedComponents::init() {
        edgesMap.clear();
        compSize.clear();
        components.assign(G->upperNodeIdBound(), none);
        tmpDistances.assign(G->upperNodeIdBound(), none);
        indexEdges();
        isTree.assign(edgesMap.size(), false);
        hasRun = false;
    }

    void DynConnectedComponents::run() {
        // Initializing / resetting data structures
        init();
        std::queue<node> q;

        // Perform breadth-first searches
        G->forNodes([&](node u) {
            if (components[u] == none) {
                index c = compSize.size();
                q.push(u);
                components[u] = c;
                compSize.insert(std::pair<index, count>(c, 1));

                do {
                    node u = q.front();
                    q.pop();
                    G->forNeighborsOf(u, [&](node v) {
                        if (components[v] == none) {
                            q.push(v);
                            components[v] = c;
                            isTree[(int)edgesMap.find(
                                makePair(u, v)
                            )->second] = true;
                            compSize.find(c)->second += 1;
                        }
                    });
                } while (!q.empty());
            }
        });

        hasRun = true;
    }


    void DynConnectedComponents::indexEdges() {
        edgeid eid = 0;
        G->forEdges([&](node u, node v) {
            if (edgesMap.find(makePair(u, v)) == edgesMap.end()) {
                insertEdgeIntoMap(u, v, eid);
                ++eid;
            }
        });
    }

    void DynConnectedComponents::insertEdgeIntoMap(node u, node v, index eid) {
        edgesMap.insert(std::pair<std::pair<node, node>, index>(
            makePair(u, v), eid)
        );
    }


    void DynConnectedComponents::update(GraphEvent event) {
        assureFinished();

        if (event.type == GraphEvent::EDGE_ADDITION) {
            addEdge(event.u, event.v);
        }
        else if (event.type == GraphEvent::EDGE_REMOVAL) {
            removeEdge(event.u, event.v);
        }
    }


    void DynConnectedComponents::updateBatch(
        const std::vector<GraphEvent>& batch
    ){
        for (auto e : batch) {
            update(e);
        }
    }


    std::pair<bool, edgeid> DynConnectedComponents::updateMapAfterAddition(
        node u, node v
    ) {
        auto it = edgesMap.find(makePair(u, v));

        if (it == edgesMap.end()) {
            edgeid newId = edgesMap.size();
            // Adding edge never deleted before
            insertEdgeIntoMap(u, v, newId);
            return std::make_pair(false, none);
        }
        return std::make_pair(true, it->second);
    }


    void DynConnectedComponents::addEdge(node u, node v) {

        std::pair<bool, edgeid> updateResult = updateMapAfterAddition(u, v);

        // If u and v are already in the same component, we
        // don't have to do anything
        index maxComp = std::max(components[u], components[v]);
        index minComp = std::min(components[u], components[v]);

        if (maxComp == minComp) {
            if (!updateResult.first) {
                isTree.push_back(false);
            }
            return;
        }

        // In the other case, we can merge the two components in an undirected
        // graph merge components
        G->parallelForNodes([&](node w) {
            if (components[w] == maxComp) {
                components[w] = minComp;
            }
        });

        compSize.find(minComp)->second += compSize.find(maxComp)->second;
        compSize.erase(maxComp);
        componentIds.push(maxComp);

        if (updateResult.first) {
            isTree[updateResult.second] = true;
        }
        else {
            isTree.push_back(true);
        }
    }


    void DynConnectedComponents::removeEdge(node u, node v) {

        edgeid eid = edgesMap.find(makePair(u, v))->second;

        // This edge removal does not split two components. Nothing to do.
        if (!isTree[eid]) {
            return;
        }

        isTree[eid] = false; // for coherence, we mark this edge as not valid
        std::fill(tmpDistances.begin(), tmpDistances.end(), none);
        index nextId = nextAvailableComponentId(false);

        std::vector<node> newCmp(components);
        newCmp[u] = nextId;
        count newCmpSize = 0;

        std::queue<node> q;
        q.push(u);
        tmpDistances[u] = 0;

        bool connected = false;

        // Berform BFS from v to check if v reaches u
        do {
            node s = q.front();
            q.pop();

            count d = tmpDistances[s] + 1;

            // Enqueue not visited neighbors
            G->forNeighborsOf(s, [&](node w) {
                if (!connected) {
                    if (tmpDistances[w] == none) {
                        tmpDistances[w] = d;
                        if (w == v) { // Found another path from u to v
                            // Backtracks the path from v to u and marks all its
                            // nodes as part of the spanning tree
                            reverseBFS(u, v);
                            connected = true;
                            return; // Exit from the loop
                        }

                        newCmp[w] = nextId;
                        ++newCmpSize;
                        q.push(w);
                    }
                }
            });

            if (connected) {
                break;
            }

        } while (!q.empty());

        if (!connected) {
            index nextId = nextAvailableComponentId();
            compSize.find(components[u])->second -= newCmpSize;
            compSize.insert(std::pair<index, count>(nextId, newCmpSize));
            components = newCmp;
        }
    }


    void DynConnectedComponents::reverseBFS(node u, node v) {

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

                    if (w == u) {
                        isTree[edgesMap.find(makePair(w, s))->second] = true;
                        nextEdgeFound = true;
                        return;
                    }

                    if (!nextEdgeFound) {
                        if ((tmpDistances[w] != none)
                        && (d == tmpDistances[w] + level)) {
                            isTree[edgesMap.find(
                                makePair(w, s)
                            )->second] = true;
                            nextEdgeFound = true;
                            q.push(w);
                        }
                    }
                }
            });
            ++level;
        } while (!q.empty());
    }


    index DynConnectedComponents::nextAvailableComponentId(bool eraseId) {
        if (componentIds.empty()) {
            return compSize.size();
        }
        index result = componentIds.front();
        if (eraseId) {
            componentIds.pop();
        }
        return result;
    }


    std::vector<std::vector<node> > DynConnectedComponents::getComponents() {
        assureFinished();

        std::vector<std::vector<node> > result(compSize.size());
        std::map<index, count> compIndex;

        int i = 0;
        for (auto it=compSize.begin(); it!=compSize.end(); ++it) {
            auto indexIterator = compIndex.find(it->first);
            if (indexIterator == compIndex.end()) {
                compIndex.insert(std::pair<index, count>(it->first, i));
                ++i;
            }
        }

        G->forNodes([&](node u) {
            result[compIndex.find(components[u])->second].push_back(u);
        });

        return result;
    }



    std::pair<node, node> DynConnectedComponents::makePair(node u, node v) {
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

}
