/*
 *
 */

#include <networkit/graph/RandomMaximumSpanningForest.hpp>
#include <networkit/auxiliary/SignalHandling.hpp>
#include <networkit/auxiliary/Parallel.hpp>
#include <networkit/graph/GraphTools.hpp>

namespace NetworKit {

RandomMaximumSpanningForest::RandomMaximumSpanningForest(const Graph &G) : G(G), hasWeightedEdges(false), hasMSF(false), hasAttribute(false) { }

void RandomMaximumSpanningForest::run() {
    hasRun = false;
    hasMSF = false;
    hasAttribute = false;

    Aux::SignalHandler handler;

    msf = GraphTools::copyNodes(G);

    handler.assureRunning();

    bool useEdgeWeights = false;

    if (!hasWeightedEdges) {
        weightedEdges.reserve(G.numberOfEdges());

        G.forEdges([&](node u, node v, edgeweight weight, edgeid eid) {
            weightedEdges.emplace_back(u, v, weight, eid);
        });

        hasWeightedEdges = true;
        useEdgeWeights = true;
    }

    handler.assureRunning();

    bool calculateAttribute = false;

    if (G.hasEdgeIds()) {
        msfAttribute.clear();
        msfAttribute.resize(G.upperEdgeIdBound(), false);
        calculateAttribute = true;
    }

    Aux::Parallel::sort(weightedEdges.begin(), weightedEdges.end(), std::greater<weightedEdge>());

    handler.assureRunning();

    UnionFind uf(G.upperNodeIdBound());

    for (weightedEdge e : weightedEdges) {
        if (uf.find(e.u) != uf.find(e.v)) {
            if (useEdgeWeights) {
                msf.addEdge(e.u, e.v, e.attribute);
            } else {
                msf.addEdge(e.u, e.v);
            }

            if (calculateAttribute) {
                msfAttribute[e.eid] = true;
            }

            uf.merge(e.u, e.v);
        }
    }

    handler.assureRunning();

    hasAttribute = calculateAttribute;
    hasMSF = true;
    hasRun = true;
}

bool RandomMaximumSpanningForest::inMSF(edgeid eid) const {
    if (!hasAttribute) throw std::runtime_error("Error: Either the attribute hasn't be calculated yet or the graph has no edge ids.");

    return msfAttribute[eid];
}

bool RandomMaximumSpanningForest::inMSF(node u, node v) const {
    if (hasMSF) {
        return msf.hasEdge(u, v);
    } else if (hasAttribute) {
        return msfAttribute[G.edgeId(u, v)];
    } else {
        throw std::runtime_error("Error: The run() method must be executed first");
    }
}

std::vector< bool > RandomMaximumSpanningForest::getAttribute(bool move) {
    std::vector<bool> result;

    if (!hasAttribute) throw std::runtime_error("Error: The run() method must be executed first");

    if (move) {
        result = std::move(msfAttribute);
        hasAttribute = false;
    } else {
        result = msfAttribute;
    }

    return result;
}

const Graph& RandomMaximumSpanningForest::getMSF() const {
    this->assureFinished();
    return msf;
}

std::string RandomMaximumSpanningForest::toString() const {
    return "Random maximum weight spanning forest";
}

bool RandomMaximumSpanningForest::isParallel() const {
    return false;
}


} // namespace NetworKit
