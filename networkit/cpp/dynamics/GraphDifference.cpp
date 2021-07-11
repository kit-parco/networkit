// no-networkit-format
#include <string>

#include <networkit/auxiliary/Log.hpp>
#include <networkit/dynamics/GraphDifference.hpp>

namespace NetworKit {

GraphDifference::GraphDifference(const Graph &G1, const Graph &G2) : G1(&G1), G2(&G2) {
    if (G1.isDirected() != G2.isDirected()) {
        throw std::runtime_error("Error, either both or none of the graphs must be directed.");
    }

    if (G1.isWeighted() != G2.isWeighted()) {
        throw std::runtime_error("Error, either both or none of the graphs must be weighted.");
    }
}

void GraphDifference::run() {
    hasRun = false;
    edits.clear();
    numNodeAdditions = 0;
    numNodeRemovals = 0;
    numNodeRestorations = 0;
    numEdgeAdditions = 0;
    numEdgeRemovals = 0;
    numWeightUpdates = 0;

    std::vector<bool> marker(G1->upperNodeIdBound(), false);
    std::vector<edgeweight> neighborWeights(G1->upperNodeIdBound(), 0);

    // collect node events and edge removals/additions in separate vectors
    // so we can later put them in the right order: first remove edges,
    // then remove and add nodes and then add edges.
    std::vector<GraphEvent> nodeEvents, edgeRemovals, edgeAdditions;

    node updatedUpperNodeIdBound = G1->upperNodeIdBound();
    for (node u = 0; u < G1->upperNodeIdBound() || u < G2->upperNodeIdBound(); ++u) {
        // First, fix non-common nodes
        if (!G2->hasNode(u) && G1->hasNode(u)) {
            nodeEvents.emplace_back(GraphEvent::NODE_REMOVAL, u);
            ++numNodeRemovals;
        } else if (G2->hasNode(u) && !G1->hasNode(u)) {
            if (u < G1->upperNodeIdBound()) {
                nodeEvents.emplace_back(GraphEvent::NODE_RESTORATION, u);
                ++numNodeRestorations;
            } else {
                while (u > updatedUpperNodeIdBound) {
                    // add and remove the same node immediately until we reach
                    // the desired node id.
                    nodeEvents.emplace_back(GraphEvent::NODE_ADDITION);
                    nodeEvents.emplace_back(GraphEvent::NODE_REMOVAL, updatedUpperNodeIdBound);
                    ++updatedUpperNodeIdBound;
                }

                // add the actually wanted node.
                nodeEvents.emplace_back(GraphEvent::NODE_ADDITION);
                ++updatedUpperNodeIdBound;
                ++numNodeAdditions;
            }
        }

        // mark neighbors of current node in G1
        if (G1->hasNode(u)) {
            G1->forNeighborsOf(u, [&](node v, edgeweight w) {
                    if (G1->isDirected() || u <= v) {
                        TRACE("Marking neighbor of ", u, " in G1 ", v);
                        marker[v] = true;
                        neighborWeights[v] = w;
                    }
                });
        }

        // unmark common neighbors, detect edge additions
        if (G2->hasNode(u)) {
            G2->forNeighborsOf(u, [&](node v, edgeweight w) {
                    // for undirected graphs, edges are only added in one direction unless
                    // the other node does not exist in G1 (edges where both nodes do not
                    // exist in G1 have been added above).
                    if (G1->isDirected() || u <= v) {
                        if (v < G1->upperNodeIdBound() && marker[v]) {
                            if (neighborWeights[v] != w) {
                                edgeAdditions.emplace_back(GraphEvent::EDGE_WEIGHT_UPDATE, u, v, w);
                                ++numWeightUpdates;
                            }
                            TRACE("Unmarking neighbor of ", u, " in G2 ", v);
                            marker[v] = false;
                        } else {
                            edgeAdditions.emplace_back(GraphEvent::EDGE_ADDITION, u, v, w);
                            ++numEdgeAdditions;
                        }
                    }
                });
        }

        if (G1->hasNode(u)) {
            // detect edge removals, unset the remaining neighbor markers
            G1->forNeighborsOf(u, [&](node v) {
                    TRACE("Checking again (", u, ",", v, ")");
                    if (G1->isDirected() || u <= v) {
                        TRACE("Edge (", u, ",", v, ") is considered");
                        if (marker[v]) {
                            TRACE("Deleting (", u, ",", v, ")");
                            edgeRemovals.emplace_back(GraphEvent::EDGE_REMOVAL, u, v);
                            ++numEdgeRemovals;
                            marker[v] = false;
                        }
                    }
                });
        }
    }

    numEdits = numNodeRemovals + numNodeAdditions + numNodeRestorations + numEdgeRemovals + numEdgeAdditions + numWeightUpdates;

    edits.swap(edgeRemovals);
    edits.reserve(nodeEvents.size() + edgeAdditions.size());
    edits.insert(edits.end(), nodeEvents.begin(), nodeEvents.end());
    edits.insert(edits.end(), edgeAdditions.begin(), edgeAdditions.end());
    hasRun = true;
}

std::vector< GraphEvent > GraphDifference::getEdits() const {
    assureFinished();

    return edits;
}

count GraphDifference::getNumberOfEdits() const {
    assureFinished();

    return numEdits;
}

count GraphDifference::getNumberOfNodeAdditions() const {
    assureFinished();

    return numNodeAdditions;
}

count GraphDifference::getNumberOfNodeRemovals() const {
    assureFinished();

    return numNodeRemovals;
}

count GraphDifference::getNumberOfNodeRestorations() const {
    assureFinished();

    return numNodeRestorations;
}

count GraphDifference::getNumberOfEdgeAdditions() const {
    assureFinished();

    return numEdgeAdditions;
}

count GraphDifference::getNumberOfEdgeRemovals() const {
    assureFinished();

    return numEdgeRemovals;
}

count GraphDifference::getNumberOfEdgeWeightUpdates() const {
    assureFinished();

    return numWeightUpdates;
}

} // namespace NetworKit
