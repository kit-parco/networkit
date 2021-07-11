/*
 * PubWebGenerator.hpp
 *
 *  Created on: Apr 10, 2013
 *      Author: Henning
 */

#ifndef NETWORKIT_GENERATORS_PUB_WEB_GENERATOR_HPP_
#define NETWORKIT_GENERATORS_PUB_WEB_GENERATOR_HPP_

#include <vector>

#include <networkit/generators/StaticGraphGenerator.hpp>
#include <networkit/graph/Graph.hpp>
#include <networkit/viz/Point.hpp>

namespace NetworKit {

/**
 * @ingroup generators
 * Generates a static graph that resembles an assumed geometric distribution of nodes in
 * a P2P network. The basic structure is to distribute points randomly in the unit torus
 * and to connect vertices close to each other (at most @a neighRad distance and none of
 * them already has @a maxNeigh neighbors). The distribution is chosen to get some areas with
 * high density and others with low density. There are @a numDenseAreas dense areas, which can
 * overlap. Each area is circular, has a certain position and radius and number of points.
 * These values are strored in @a denseAreaXYR and @a numPerArea, respectively.
 *
 * Used and described in more detail in J. Gehweiler, H. Meyerhenke: A Distributed
 * Diffusive Heuristic for Clustering a Virtual P2P Supercomputer. In Proc. 7th High-Performance
 * Grid Computing Workshop (HPGC'10), in conjunction with 24th IEEE Internatl. Parallel and
 * Distributed Processing Symposium (IPDPS'10), IEEE, 2010.
 *
 * Reasonable parameters for constructor:
 * - numNodes: from several hundred up to a few thousand
 *   (possibly more if visualization is not desired and quadratic time
 *   complexity has been resolved)
 * - numberOfDenseAreas: depending on number of nodes, e.g. [8, 50]
 * - neighborhoodRadius: the higher, the better the connectivity [0.1, 0.35]
 * - maxNumberOfNeighbors: maximum degree, a higher value corresponds to better connectivity [4, 40]
 */
class PubWebGenerator final : public StaticGraphGenerator {
    friend class DynamicPubWebGenerator;

public:
    PubWebGenerator() {
    } // nullary constructor needed for Python Shell - do not use this to construct instance
    ~PubWebGenerator() override = default;

    PubWebGenerator(count numNodes, count numberOfDenseAreas, coordinate neighborhoodRadius,
                    count maxNumberOfNeighbors);

    Graph generate() override;

    const std::vector<Point2D> &getCoordinates() const { return coordinates; }
    std::vector<Point2D> moveCoordinates() { return std::move(coordinates); }

private:
    struct circle {
        coordinate x;
        coordinate y;
        coordinate rad;
    };

    static constexpr coordinate MAX_DENSE_AREA_RADIUS = 0.2;
    static constexpr coordinate MIN_MAX_DENSE_AREA_FACTOR = 5.0;
    static constexpr edgeweight BASE_WEIGHT = 0.01;

    count n;                          //!< number of nodes
    count numDenseAreas;              //!< number of areas with more nodes (denser)
    coordinate neighRad;              //!< neighborhood radius
    count maxNeigh;                   //!< maximum number of neighbors
    std::vector<circle> denseAreaXYR; //!< position of each circular dense area
    std::vector<count> numPerArea;    //!< number of points in each circular area
    std::vector<Point2D> coordinates; //!< storage for point coordinates

    Point2D intoUnitSquare(Point2D pt) const noexcept;
    coordinate squaredDistanceInUnitTorus(Point2D pt1, Point2D pt2) const noexcept;

    void determineNeighbors(Graph &g);
    void chooseDenseAreaSizes();
    void fillDenseAreas(Graph &g);
    void spreadRemainingNodes(Graph &g);
    void chooseClusterSizes();
    void addNodesToArea(index area, count num, Graph &g);

    /**
     * Adds nodes randomly, distribution respects original one.
     */
    void addNode(Graph &g);
};

} /* namespace NetworKit */
#endif // NETWORKIT_GENERATORS_PUB_WEB_GENERATOR_HPP_
