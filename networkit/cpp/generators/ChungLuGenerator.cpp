// no-networkit-format
/*
 * ChungLu.cpp
 *
 *  Created on: Dec 23, 2013
 *      Author: Henning
 *      Contributors: Hoske/Weisbarth
 */

#include <numeric>

#include <networkit/auxiliary/Parallel.hpp>
#include <networkit/auxiliary/Random.hpp>
#include <networkit/generators/ChungLuGenerator.hpp>
#include <networkit/graph/GraphBuilder.hpp>

namespace NetworKit {

ChungLuGenerator::ChungLuGenerator(const std::vector<count> &degreeSequence) :
        StaticDegreeSequenceGenerator(degreeSequence) {
    sum_deg = std::accumulate(seq.begin(), seq.end(), count{0});
    n = (count) seq.size();
}

    Graph ChungLuGenerator::generate() {
        GraphBuilder gB(n);

        /* We need a sorted list in descending order for this algorithm */
        Aux::Parallel::sort(seq.begin(), seq.end(), [](count a, count b){ return a > b;});

        for (node u = 0; u <= n - 2; u++) {
            node v = u + 1;
            /* Apparently it is necessary to include all these casts for
             * the probability to be properly calculated */
            double p = std::min(((double) seq[u]) * ((double) seq[v]) / sum_deg, 1.0);

            while (v < n && p > 0) {
                if (p != 1.0) {
                    double randVal = Aux::Random::probability();
                    /* Calculate the distance to the next potential neighbour*/
                    v = v + (node) std::floor(log(randVal)/log(1 - p));
                }
                if ((count) v < n) {
                    double q = std::min(((double) seq[u]) * ((double) seq[v]) / sum_deg, 1.0);
                    double randVal2 = Aux::Random::probability();
                    /* The potential neighbor was selected with the probability p.
                     * In order to see if this neighbor should be rejected or accepted
                     * we correct the probability using q */
                    if (randVal2 < q / p) {
                        gB.addHalfOutEdge(u, v);
                    }
                    p = q;
                    v++;
                }
            }
        }

        return gB.toGraph(true,true);
    }

} /* namespace NetworKit */
