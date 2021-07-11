// no-networkit-format
/*
 * EdgeScoreLinearizer.cpp
 *
 *  Created on: 18.11.2014
 *      Author: Michael Hamann
 */

#include <numeric>
#include <tuple>

#include <networkit/auxiliary/Parallel.hpp>
#include <networkit/auxiliary/Random.hpp>
#include <networkit/edgescores/EdgeScoreLinearizer.hpp>

namespace NetworKit {

EdgeScoreLinearizer::EdgeScoreLinearizer(const Graph& G, const std::vector<double>& attribute, bool inverse) : EdgeScore<double>(G), attribute(&attribute), inverse(inverse) {
}

void EdgeScoreLinearizer::run() {
    if (!G->hasEdgeIds()) {
        throw std::runtime_error("edges have not been indexed - call indexEdges first");
    }

    scoreData.resize(G->upperEdgeIdBound());

    // Special case for m = 1
    if (G->numberOfEdges() == 1) {
        G->forEdges([&](node, node, edgeid eid) {
            scoreData[eid] = 0.5;
        });
    } else {
        using edgeTuple = std::tuple<edgeweight, index, edgeid>;
        std::vector<edgeTuple> sorted(G->upperEdgeIdBound(), std::make_tuple(std::numeric_limits<edgeweight>::max(), std::numeric_limits<index>::max(), none));

        G->parallelForEdges([&](node, node, edgeid eid) {
            sorted[eid] = std::make_tuple((*attribute)[eid], Aux::Random::integer(), eid);
        });

        if (inverse) {
            Aux::Parallel::sort(sorted.begin(), sorted.end(), std::greater<edgeTuple>());
        } else {
            Aux::Parallel::sort(sorted.begin(), sorted.end());
        }

        #pragma omp parallel for
        for (omp_index pos = 0; pos < static_cast<omp_index>(G->upperEdgeIdBound()); ++pos) {
            edgeid eid = std::get<2>(sorted[pos]);
            if (eid != none) {
                scoreData[eid] = static_cast<double>(pos) * 1.0 / static_cast<double>(G->numberOfEdges()-1); // TODO maybe writing back to sorted and then sorting again by eid could be faster (cache efficiency)?
            }
        }
    }

    hasRun = true;
}

double EdgeScoreLinearizer::score(node, node) {
    throw std::runtime_error("Not implemented: Use scores() instead.");
}

double EdgeScoreLinearizer::score(edgeid) {
    throw std::runtime_error("Not implemented: Use scores() instead.");
}

} /* namespace NetworKit */
