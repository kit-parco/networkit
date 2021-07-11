// no-networkit-format
/*
 * PivotMDS.hpp
 *
 *  Created on: Jul 7, 2016
 *      Author: Michael Wegner
 */

#ifndef NETWORKIT_VIZ_PIVOT_MDS_HPP_
#define NETWORKIT_VIZ_PIVOT_MDS_HPP_

#include <networkit/algebraic/CSRMatrix.hpp>
#include <networkit/algebraic/Vector.hpp>
#include <networkit/graph/Graph.hpp>
#include <networkit/viz/GraphLayoutAlgorithm.hpp>

namespace NetworKit {

/**
 * @ingroup viz
 *
 * Implementation of PivotMDS proposed by Brandes and Pich.
 */
class PivotMDS final : public GraphLayoutAlgorithm<double> {
  public:
    /**
     * Constructs a PivotMDS object for the given @a graph. The algorithm should
     * embed the graph in @a dim dimensions using @a numPivots pivots.
     * @param graph
     * @param dim
     * @param numPivots
     */
    PivotMDS(const Graph &graph, count dim, count numPivots);

    /*
     * Default destructor
     */
    ~PivotMDS() override = default;

    /**
     * Runs the PivotMDS algorithm.
     */
    void run() override;

  private:
    count dim;
    count numPivots;

    /**
     *  Randomly picks the pivots for the algorithm
     */
    std::vector<node> computePivots();

    /**
     * Power method to compute the largest eigenvector and eigenvalue that are
     * stored in @a eigenvector and
     * @a eigenvalue.
     */
    void powerMethod(const CSRMatrix &mat, count n, Vector &eigenvector,
                     double &eigenvalue);
};

} /* namespace NetworKit */

#endif // NETWORKIT_VIZ_PIVOT_MDS_HPP_
