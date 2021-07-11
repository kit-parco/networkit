// no-networkit-format
/*
 * PrefixJaccardScore.cpp
 *
 *  Created on: 26.09.2014
 *      Author: Michael Hamann
 */

#include <omp.h>

#include <networkit/edgescores/PrefixJaccardScore.hpp>

namespace NetworKit {

template <typename AttributeT>
PrefixJaccardScore<AttributeT>::PrefixJaccardScore(const Graph &G, const std::vector< AttributeT > &attribute) :
    EdgeScore<double>(G), inAttribute(&attribute) {
}

template <typename AttributeT>
void PrefixJaccardScore<AttributeT>::run() {
    //this-> required to access members of base class, since this is a template class.
    if (!this->G->hasEdgeIds()) throw std::runtime_error("Error, edges need to be indexed first");

    this->scoreData.clear();
    this->scoreData.resize(this->G->upperEdgeIdBound());

    struct RankedEdge {
        node u;
        AttributeT att;
        count rank;

        RankedEdge(node u, AttributeT att, count rank) : u(u), att(att), rank(rank) {};

        bool operator<(const RankedEdge &other) const {
            return std::tie(rank, att, u) < std::tie(other.rank, other.att, other.u);
        };

        bool operator>(const RankedEdge &other) const {
            return std::tie(rank, att, u) > std::tie(other.rank, other.att, other.u);
        };
    };

    std::vector<size_t> rankedEdgeBegin(G->upperNodeIdBound() + 1);
    std::vector<RankedEdge> rankedEdges;
    rankedEdges.reserve(2*G->numberOfEdges());

    for (node u = 0; u < G->upperNodeIdBound(); ++u) {
        rankedEdgeBegin[u] = rankedEdges.size();
        if (G->hasNode(u)) {
            G->forEdgesOf(u, [&](node, node v, edgeid eid) {
                rankedEdges.emplace_back(v, ((*inAttribute)[eid]), 0);
            });
        }
    }
    rankedEdgeBegin[G->upperNodeIdBound()] = rankedEdges.size();

    this->G->balancedParallelForNodes([&](node u) {
        if (this->G->degree(u) == 0) return;

        const auto beginIt = rankedEdges.begin() + rankedEdgeBegin[u];
        const auto endIt = rankedEdges.begin() + rankedEdgeBegin[u+1];

        std::sort(beginIt, endIt, std::greater<RankedEdge>());

        AttributeT curVal = beginIt->att;
        count curRank = 0;
        count numEqual = 0;
        for (auto it = beginIt; it != endIt; ++it) {
            if (curVal != it->att) {
                curRank += numEqual;
                curVal = it->att;
                numEqual = 1;
            } else {
                ++numEqual;
            }

            it->rank = curRank;
        }
    });

    std::vector<std::vector<bool>> uMarker(omp_get_max_threads(), std::vector<bool>(G->upperNodeIdBound(), false));
    auto vMarker = uMarker;

    this->G->parallelForEdges([&](node u, node v, edgeid eid) {
        count curRank = 0;
        double bestJaccard = 0;
        auto tid = omp_get_thread_num();

        auto uIt = rankedEdges.begin() + rankedEdgeBegin[u];
        auto vIt = rankedEdges.begin() + rankedEdgeBegin[v];
        const auto uEndIt = rankedEdges.begin() + rankedEdgeBegin[u+1];
        const auto vEndIt = rankedEdges.begin() + rankedEdgeBegin[v+1];

        double commonNeighbors = 0.0;
        double uNeighbors = 0.0;
        double vNeighbors = 0.0;

        while (uIt != uEndIt || vIt != vEndIt) {
            while (uIt != uEndIt && curRank == uIt->rank) {
                if (uIt->u == v) {
                    ++uIt;
                    continue;
                }

                if (vMarker[tid][uIt->u]) {
                    vMarker[tid][uIt->u] = false;
                    ++commonNeighbors;
                    --vNeighbors;
                } else {
                    uMarker[tid][uIt->u] = true;
                    ++uNeighbors;
                }

                ++uIt;
            }

            while (vIt != vEndIt && curRank == vIt->rank) {
                if (vIt->u == u) {
                    ++vIt;
                    continue;
                }

                if (uMarker[tid][vIt->u]) {
                    uMarker[tid][vIt->u] = false;
                    ++commonNeighbors;
                    --uNeighbors;
                } else {
                    vMarker[tid][vIt->u] = true;
                    ++vNeighbors;
                }

                ++vIt;
            }

            bestJaccard = std::max(bestJaccard, commonNeighbors * 1.0 / (uNeighbors + vNeighbors + commonNeighbors));

            ++curRank;
        }

        G->forNeighborsOf(u, [&](node w) {
            uMarker[tid][w] = false;
        });

        G->forNeighborsOf(v, [&](node w) {
            vMarker[tid][w] = false;
        });

        this->scoreData[eid] = bestJaccard;
    });

    this->hasRun = true;
}

template <typename AttributeT>
double PrefixJaccardScore<AttributeT>::score(node, node) {
    throw std::runtime_error("Not implemented: Use scores() instead.");
}

template <typename AttributeT>
double PrefixJaccardScore<AttributeT>::score(edgeid) {
    throw std::runtime_error("Not implemented: Use scores() instead.");
}

template class PrefixJaccardScore<double>;
template class PrefixJaccardScore<count>;

}
