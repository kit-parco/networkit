// no-networkit-format
/*
 * SimmelianScore.cpp
 *
 *  Created on: 21.05.2014
 *      Author: Gerd Lindner
 */

#include <networkit/sparsification/SimmelianScore.hpp>
#include <limits>

namespace NetworKit {

SimmelianScore::SimmelianScore(const Graph& G, const std::vector<count>& triangles) : EdgeScore<double>(G), triangles(&triangles) {
}

std::vector<RankedNeighbors> SimmelianScore::getRankedNeighborhood(const Graph& g, const std::vector<count>& triangles) {
    std::vector<RankedNeighbors> neighbors;
    neighbors.resize(g.upperNodeIdBound());

    g.forNodes([&](node u) {
        //Sort ego's alters from strongly to weakly tied.
        g.forNeighborsOf(u, [&](node, node v, edgeid eid) {
            count triangleCount = round(triangles[eid]);
            neighbors[u].push_back(RankedEdge(u, v, triangleCount));
        });
        std::sort(neighbors[u].begin(), neighbors[u].end());

        //Calculate the ranks.
        count currentRank = 0;	//Rank 0 is considered the best.
        count currentSimmelianness = std::numeric_limits<count>::max();
        count equals = 0;
        for (auto& edge : neighbors[u]) {
            if (edge.simmelianness != currentSimmelianness) {
                currentRank += equals;
                currentSimmelianness = edge.simmelianness;
                equals = 1;
            } else {
                equals++;
            }
            edge.rank = currentRank;
        }
    });

    return neighbors;

}

Redundancy SimmelianScore::getOverlap(	const node& ego,
                                                const node& alter,
                                                const std::vector<RankedNeighbors>& neighbors,
                                                const count& maxRank) {
    //Initialization of output values
    Redundancy result = Redundancy(0, 0.0);

    std::vector<RankedEdge>::const_iterator egoIt = neighbors[ego].begin();
    std::vector<RankedEdge>::const_iterator alterIt = neighbors[alter].begin();

    std::set<node> egoNeighborsUnmatched;
    std::set<node> alterNeighborsUnmatched;

    for (count rank = 0; rank <= maxRank; rank++) {
        matchNeighbors(ego, alter, true, egoIt, neighbors[ego], egoNeighborsUnmatched, alterNeighborsUnmatched, rank, result.overlap);
        matchNeighbors(alter, ego, false, alterIt, neighbors[alter], alterNeighborsUnmatched, egoNeighborsUnmatched, rank, result.overlap);

        double currentJaccard = 0.0;
        if (result.overlap + egoNeighborsUnmatched.size() + alterNeighborsUnmatched.size() > 0)
            currentJaccard = double(result.overlap) / double(result.overlap + egoNeighborsUnmatched.size() + alterNeighborsUnmatched.size());

        result.jaccard = std::max(currentJaccard, result.jaccard);
    }

    return result;
}

/**
 * Helper function used in getOverlap. Adds the intersection of
 * egoNeighbors and alterNeighborsUnmatched to overlap.
 */
void SimmelianScore::matchNeighbors(node, node alter, bool,
                                    std::vector<RankedEdge>::const_iterator &,
                                    const RankedNeighbors &egoNeighbors,
                                    std::set<node> &egoNeighborsUnmatched,
                                    std::set<node> &alterNeighborsUnmatched, count rank,
                                    count &overlap) {

    for (auto egoIt : egoNeighbors) {
        node other = egoIt.alter;

        if (other == alter || egoIt.rank != rank)
            continue;

        if (alterNeighborsUnmatched.erase(other))
            overlap++;
        else
            egoNeighborsUnmatched.insert(other);
    }
}

double SimmelianScore::score(node, node) {
    throw std::runtime_error("Not implemented: Use scores() instead.");
}

double SimmelianScore::score(edgeid) {
    throw std::runtime_error("Not implemented: Use scores() instead.");
}

} /* namespace NetworKit */
