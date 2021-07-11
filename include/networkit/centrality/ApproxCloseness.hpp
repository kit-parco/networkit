// no-networkit-format
/*
 * ApproxCloseness.h
 *
 *  Created on: Dec 8, 2015
 *      Author: Sarah Lutteropp (uwcwa@student.kit.edu) and Michael Wegner (michael.wegner@student.kit.edu)
 */

#ifndef NETWORKIT_CENTRALITY_APPROX_CLOSENESS_HPP_
#define NETWORKIT_CENTRALITY_APPROX_CLOSENESS_HPP_

#include <networkit/centrality/Centrality.hpp>
#include <limits>

namespace NetworKit {

/**
 * @ingroup centrality
 * Approximation of closeness centrality according to algorithm described in
 * Cohen et al., Computing Classic Closeness Centrality, at Scale
 */
class ApproxCloseness: public Centrality {

public:
    enum CLOSENESS_TYPE {INBOUND, OUTBOUND, SUM};

    /**
     * Constructs an instance of the ApproxCloseness class for @a graph using @a nSamples during the run() method.
     * The @a epsilon parameter (standard = 0.1) is used to control the switch between sampling and pivoting.
     * Using @a epsilon = 0, the algorithm only uses sampling. (see Cohen, Edith, et al.
     * "Computing classic closeness centrality, at scale." Proceedings of the second ACM conference on Online social
     * networks. ACM, 2014.). The running time is proportional to nSamples * m, where m is  the number of edges.
     * Notice: the input graph has to be connected.
     * @param	graph		input graph
     * @param	nSamples	user defined number of samples
     * @param 	epsilon		Value in [0, infty) controlling the switch between sampling and pivoting. When using 0, only sampling is used. Standard is 0.1.
     * @param	normalized  normalize centrality values in interval [0,1]
     * @param 	type		use in- or outbound centrality or the sum of both (see paper) for computing closeness on directed graph. If G is undirected, this can be ignored.
     */
    ApproxCloseness(const Graph& G, count nSamples, double epsilon = 0.1, bool normalized=false, CLOSENESS_TYPE type = OUTBOUND);


    /**
     * Computes closeness approximation on the graph passed in constructor.
     */
    void run() override;

    /**
     * Returns the maximum possible Closeness a node can have in a graph with the same amount of nodes (=a star)
     */
    double maximum() override;

    /**
     * @return The square error when closeness centrality has been computed for an undirected graph.
     */
    std::vector<double> getSquareErrorEstimates();

private:
    count nSamples;
    double epsilon;

    /** \sum_{i \in L(j) \cap C} d_{ji} for every j \in V */
    std::vector<double> LCSum;

    /** |L(j) \cap C| for every j \in V */
    std::vector<count> LCNum;

    /** \sum_{i \in L(j) \cap C} d_{ji}^2 for every j \in V */
    std::vector<double> LCSumSQ;

    /** \sum_{i \in HC(j)} d_{ji} for every j \in V */
    std::vector<double> HCSum;

    /** \frac{1}{|HC(j)|} * \sum_{i \in HC(j)} (d_{ij} - d_{c(j)j})^2 for every j \in V */
    std::vector<double> HCSumSQErr;

    /** \sum_{i \in H(j)} d_{c(j)i} for every j \in V */
    std::vector<double> HSum;

    /** |H(j)| for every j \in V */
    std::vector<count> HNum;

    /** Reachability estimation **/
    std::vector<double> R;

    std::vector<double> SQErrEst;
    const edgeweight infDist = floor(std::numeric_limits<edgeweight>::max() / 2.0); // divided by two s.t. infDist + infDist produces no overflow

    CLOSENESS_TYPE type;

    void estimateClosenessForUndirectedGraph();
    void estimateClosenessForDirectedGraph(bool outbound);
    inline void computeClosenessForDirectedWeightedGraph(bool outbound);
    inline void computeClosenessForDirectedUnweightedGraph(bool outbound);

    /**
     * Runs a multi-source Dijkstra from all @a samples and sets the closest pivot and the distance from it for
     * every node in G.
     * @param samples The sampled nodes that are the pivot nodes.
     * @param pivot[out] Stores the closest pivot (sample) from every node.
     * @param delta[out] Stores the distance d(u,pivot(u)) for every node u.
     */
    void computeClosestPivot(const std::vector<node> &samples, std::vector<node> &pivot, std::vector<edgeweight> &delta);

    /**
     * Runs the algorithm of Cohen et al. for a single pivot.
     * @param i
     * @param pivot
     * @param delta
     * @param samples
     */
    void runOnPivot(index i, const std::vector<node> &pivot, const std::vector<edgeweight> &delta, const std::vector<node> &samples);

    /**
     * Orders the nodes of G in increasing distance from @a pivot. The distances from @a pivot are stored in @a pivotDist
     * @param pivot
     * @param order[out]
     * @param pivotDist[out]
     */
    void orderNodesByIncreasingDistance(node pivot, std::vector<node> &order, std::vector<edgeweight> &pivotDist);

};

} /* namespace NetworKit */

#endif // NETWORKIT_CENTRALITY_APPROX_CLOSENESS_HPP_
