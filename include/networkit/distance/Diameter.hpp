// no-networkit-format
/*
 * Diameter.hpp
 *
 *  Created on: 19.02.2014
 *      Author: Daniel Hoske, Christian Staudt
 */

#ifndef NETWORKIT_DISTANCE_DIAMETER_HPP_
#define NETWORKIT_DISTANCE_DIAMETER_HPP_

#include <networkit/auxiliary/SignalHandling.hpp>
#include <networkit/base/Algorithm.hpp>
#include <networkit/graph/Graph.hpp>

namespace NetworKit {
enum DiameterAlgo {automatic = 0, exact = 1, estimatedRange = 2, estimatedSamples = 3, estimatedPedantic = 4};

/**
 * @ingroup distance
 */
class Diameter final : public Algorithm {

public:
    Diameter(const Graph& G, DiameterAlgo algo = DiameterAlgo::automatic, double error = -1.f, count nSamples = 0);

    void run() override;

    std::string TLX_DEPRECATED(toString() const override);

    std::pair<count, count> getDiameter() const;


private:
    const Graph* G;
    DiameterAlgo algo;
    double error;
    count nSamples;
    std::pair<count, count> diameterBounds;

    /**
     * Get the estimation of the diameter of the graph @a G. The algorithm is based on the ExactSumSweep algorithm presented in
     * Michele Borassi, Pierluigi Crescenzi, Michel Habib, Walter A. Kosters, Andrea Marino, Frank W. Takes,
     * Fast diameter and radius BFS-based computation in (weakly connected) real-world graphs: With an application to the six degrees of separation games,
     * Theoretical Computer Science, Volume 586, 27 June 2015, Pages 59-80, ISSN 0304-3975,
     * http://dx.doi.org/10.1016/j.tcs.2015.02.033.
     * (http://www.sciencedirect.com/science/article/pii/S0304397515001644)
     * @param G The graph.
     * @param error The maximum allowed relative error. Set to 0 for the exact diameter.
     * @return Pair of lower and upper bound for diameter.
     */
    std::pair<edgeweight, edgeweight> estimatedDiameterRange(const Graph& G, double error);

    /**
     * Get the exact diameter of the graph @a G. The algorithm for unweighted graphs is the same as
     * the algorithm for the estimated diameter range with error 0.
     *
     * @param G The graph.
     * @return exact diameter of the graph @a G
     */
    edgeweight exactDiameter(const Graph& G);

    /**
     * Get a 2-approximation of the node diameter (unweighted diameter) of @a G.
     *
     * @param[in] G        The graph.
     * @param[in] samples  One sample is enough if the graph is connected. If there
     *       are multiple connected components, then the number of samples
     *       must be chosen so that the probability of sampling the component
     *       with the largest diameter is high.
     * @return A 2-approximation of the vertex diameter (unweighted diameter) of @a G.
     */
    edgeweight estimatedVertexDiameter(const Graph& G, count samples);

    /** @return a 2-approximation of the vertex diameter (unweighted diameter) of @a G.
            Considers each connected component and returns the maximum diameter.
     */
    edgeweight estimatedVertexDiameterPedantic(const Graph& G);

    /** @return a 2-approximation of the vertex diameter (unweighted diameter) of @a G.
            Considers each connected component and returns the maximum diameter.
    */
    edgeweight estimatedVertexDiameterPedantic2(const Graph& G);
};

} /* namespace NetworKit */

#endif // NETWORKIT_DISTANCE_DIAMETER_HPP_
