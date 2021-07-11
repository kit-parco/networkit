// no-networkit-format
/*
 * CoreDecomposition.hpp
 *
 *  Created on: Oct 28, 2013
 *      Author: Lukas Barth, David Weiss, Christian Staudt
 */

#ifndef NETWORKIT_CENTRALITY_CORE_DECOMPOSITION_HPP_
#define NETWORKIT_CENTRALITY_CORE_DECOMPOSITION_HPP_

#include <vector>
#include <string>
#include <list>

#include <networkit/graph/Graph.hpp>
#include <networkit/centrality/Centrality.hpp>
#include <networkit/structures/Partition.hpp>
#include <networkit/structures/Cover.hpp>

#include <tlx/define/deprecated.hpp>

namespace NetworKit {

/**
 * @ingroup centrality
 * Computes k-core decomposition of a graph.
 */
class CoreDecomposition : public Centrality  {

public:

    /**
     * Create CoreDecomposition class for graph @a G. The graph may not contain self-loops.
     *
     * Contains the parallel algorithm by
     * Dasari, N.S.; Desh, R.; Zubair, M., "ParK: An efficient algorithm for k-core decomposition on multicore processors," in Big Data (Big Data), * 2014 IEEE International Conference.
     *
     * TODO complexity?
     * @param G The graph.
     * @param normalized If set to @c true the scores are normalized in the interval [0,1].
     * @param enforceBucketQueueAlgorithm If set to @c true, uses a bucket priority queue data structure. This it is generally slower than ParK but may be more flexible. TODO check
     * @param storeNodeOrder If set to @c true, the order of the nodes in ascending order of the cores is stored and can later be returned using getNodeOrder(). Enforces the sequential bucket priority queue algorithm.
     *
     * The algorithm runs in parallel if the usage of a bucket priority queue is not enforced and
     * if the node ids of the input graph are continuous (i.e., numberOfNodes() = upperNodeIdBound()).
     */
    CoreDecomposition(const Graph& G, bool normalized=false, bool enforceBucketQueueAlgorithm = false, bool storeNodeOrder = false);

    /**
     * Perform k-core decomposition of graph passed in constructor.
     */
    void run() override;

    /**
     * Get the k-cores as a graph cover object.
     *
     * @return the k-cores as a Cover
     */
    Cover getCover() const;

    /**
     * Get the k-shells as a partition object
     *
     * @return the k-shells as a Partition
     */
    Partition getPartition() const;

    /**
     * Get maximum core number.
     *
     * @return The maximum core number
     */
    index maxCoreNumber() const;

    /**
    * Get the theoretical maximum of centrality score in the given graph.
    *
    * @return The theoretical maximum centrality score.
    */
    double maximum() override;

    /**
     * Get the node order.
     *
     * This is only possible when storeNodeOrder was set.
     *
     * @return The nodes sorted by increasing core number.
     */
    const std::vector<node>& getNodeOrder() const;

    /**
     * The algorithm ParK can run in parallel under certain conditions,
     * the bucket PQ based one cannot.
     */
    bool TLX_DEPRECATED(isParallel() const override) { return canRunInParallel; }

private:

    index maxCore; // maximum core number of any node in the graph

    bool enforceBucketQueueAlgorithm; // in case one wants to switch to the alternative algorithm

    bool canRunInParallel; // signifies if a parallel algorithm can be used

    bool storeNodeOrder; // signifies if the node order shall be stored

    std::vector<node> nodeOrder; // Stores the node order, i.e., all nodes sorted by core number

    /**
     * Perform k-core decomposition of graph passed in constructor.
     * ParK is an algorithm by Naga Shailaja Dasari, Ranjan Desh, and Zubair M.
     * See http://ieeexplore.ieee.org/stamp/stamp.jsp?tp=&arnumber=7004366 for details.
     */
    void runWithParK();

    /**
     * Perform k-core decomposition of graph passed in constructor.
     * The algorithm is based on a bucket priority queue data structure.
     * It is generally slower than ParK but may be more flexible.
     */
    void runWithBucketQueues();

    /**
     * Determines nodes whose remaining degree equals @a level.
     * @param[in] level Shell number (= level) currently processed.
     * @param[in] degrees Remaining degree for each node.
     * @param[inout] curr Nodes to be processed in current level.
     */
    void scan(index level, const std::vector<count>& degrees, std::vector<node>& curr);

    /**
     * Determines in parallel the nodes whose remaining degree equals @a level.
     * @param[in] level Shell number (= level) currently processed.
     * @param[in] degrees Remaining degree for each node.
     * @param[inout] curr Nodes to be processed in current level.
     */
    void scanParallel(index level, const std::vector<count>& degrees, std::vector<node>& curr, std::vector<char>& active);

    /**
     * Processes nodes (and their neighbors) identified by previous scan.
     * @param[in] level Shell number (= level) currently processed.
     * @param[inout] degrees Remaining degree for each node.
     * @param[in] curr Nodes to be processed in this call.
     * @param[inout] next Nodes to be processed next in current level (certain neighbors of nodes in curr).
     */
    void processSublevel(index level, std::vector<count>& degrees, const std::vector<node>& curr, std::vector<node>& next);

    /**
     * Processes in parallel nodes (and their neighbors) identified by previous scan.
     * @param[in] level Shell number (= level) currently processed.
     * @param[inout] degrees Remaining degree for each node.
     * @param[in] curr Nodes to be processed in this call.
     * @param[inout] next Nodes to be processed next in current level (certain neighbors of nodes in curr).
     */
    void processSublevelParallel(index level, std::vector<count>& degrees, const std::vector<node>& curr, std::vector<node>& next, std::vector<char>& active);
};

} /* namespace NetworKit */
#endif // NETWORKIT_CENTRALITY_CORE_DECOMPOSITION_HPP_
