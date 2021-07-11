// no-networkit-format
/*
 * MaxentStress.hpp
 *
 *  Created on: 22.01.2014
 *      Author: Henning Meyerhenke and Michael Wegner
 */

#ifndef NETWORKIT_VIZ_MAXENT_STRESS_HPP_
#define NETWORKIT_VIZ_MAXENT_STRESS_HPP_

#include <memory>

#include <networkit/algebraic/CSRMatrix.hpp>
#include <networkit/auxiliary/StringTools.hpp>
#include <networkit/distance/AlgebraicDistance.hpp>
#include <networkit/numerics/LinearSolver.hpp>
#include <networkit/viz/GraphLayoutAlgorithm.hpp>
#include <networkit/viz/Octree.hpp>


namespace NetworKit {

    using CoordinateVector = std::vector<Vector>; // more meaningful and shorter name for coordinates stored by dimension

    struct ForwardEdge {
        node head;
        edgeweight weight;
    };


    /**
     * @ingroup viz
     *
     * Implementation of MaxentStress by Ganser et al. using a laplacian system solver.
     *
     * @see Ganser, Emden R., Yifan Hu and Steve North. "A maxentstress model for graph layout." Visualisation and Computer Graphics, IEEE Transsactions on 19, no. 6 (2013): 927-940.
     */

class MaxentStress final : public GraphLayoutAlgorithm<double> {
    public:
        enum GraphDistance {
            EDGE_WEIGHT,
            ALGEBRAIC_DISTANCE
        };

        enum LinearSolverType {
            LAMG,
            CONJUGATE_GRADIENT_IDENTITY_PRECONDITIONER,
            CONJUGATE_GRADIENT_DIAGONAL_PRECONDITIONER
        };

        struct ResultStats {
            double rhsTime = 0.0;
            double approxEntropyTerm = 0.0;
            double solveTime = 0.0;
        };

        MaxentStress(const Graph &G, count dim, count k, double tolerance,
                     LinearSolverType linearSolverType = LAMG, bool fastComputation = false,
                     GraphDistance graphDistance = EDGE_WEIGHT);

        MaxentStress(const Graph &G, count dim, const std::vector<Point<double>> &coordinates,
                     count k, double tolerance, LinearSolverType linearSolverType = LAMG,
                     bool fastComputation = false, GraphDistance graphDistance = EDGE_WEIGHT);

        /** Default destructor. */
        ~MaxentStress() override = default;

        /** Computes a graph drawing according to the Maxent-Stress model. */
         void run() override;

        /**
         * Scale the layout computed by @ref run() by a scalar s to minimize \sum_{u,v \in V} w_{uv} (s ||x_u - x_v|| - d_{uv}||)^2
         */
        void scaleLayout();

        /**
         * Computes a scalar s s.t. \sum_{u,v \in V} w_{uv} (s ||x_u - x_v|| - d_{uv}||)^2 is minimized.
         */
        double computeScalingFactor();

        /**
         * Computes the full stress measure of the computed layout with @ref run(). Code taken and adapted from KaDraw by Christian Schulz.
         * @note KaDraw: http://algo2.iti.kit.edu/kadraw/
         */
        double fullStressMeasure();

        /**
         * Computes the maxent stress measure for the computed layout with @a run(). Code taken and adapted from KaDraw by Christian Schulz.
         * @note KaDraw: http://algo2.iti.kit.edu/kadraw/
         */
        double maxentMeasure();

        double meanDistanceError();

        double ldme();

        /**
         * Set parameter @a q.
         * @param q
         */
        void setQ(double q) {
            this->q = q;
        }

        /**
         * Set parameter @a alpha.
         * @param alpha
         */
        void setAlpha(double alpha) {
            this->alpha = alpha;
        }

        /**
         * Set parameter @a alphaReduction.
         * @param alphaReduction
         */
        void setAlphaReduction(double alphaReduction) {
            this->alphaReduction = alphaReduction;
        }

        /**
         * Set parameter @a finalAlpha.
         * @param finalAlpha
         */
        void setFinalAlpha(double finalAlpha) {
            this->finalAlpha = finalAlpha;
        }

        void setConvergenceThreshold(double convThreshold) {
            this->convThreshold = convThreshold * convThreshold;
        }

        double getRhs() {
             return this->resultStats.rhsTime;
        };

        double getApproxEntropyTerm() {
             return this->resultStats.approxEntropyTerm;
        };

        double getSolveTime() {
            return this->resultStats.solveTime;
        };

    private:

     /**
      * Calls run() and returns the corresponding stats.
      */

     ResultStats runAlgo();

        /**
         * Reference to the linear solver to use during the maxent-stress algorithm.
         */
        LinearSolver<CSRMatrix>& solver;

        /**
         * Stats of the last time run() was called.
         */

        ResultStats resultStats;

        /** Parameters of the MaxentStress model **/
        double q, alpha, alphaReduction, finalAlpha, convThreshold;

        /** Specifies whether initial coordinates have been provided in the constructor */
        bool coordinatesProvided;

        /** Defines whether the algorithm stops when converged on a higher than the lowest level. This saves some time but
         *  usually leads to slightly worse results.
         */
        bool fastComputation;

        /** Maximum number of solves for the same value of alpha **/
        count maxSolvesPerAlpha;

        /** set known distances (S in the paper of Gansner et al.) */
        std::vector<std::vector<ForwardEdge>> knownDistances;

        /** cardinality of S */
        count knownDistancesCardinality;

        /** points of vertices are in R^{dim} */
        count dim;

        /**
         * Determines whether the run() method has already been called.
         */
        bool hasRun;

        /**
         * Checks whether the MaxentStress algorithm converged, i.e. ||newCoords - oldCoords|| / ||oldCoords|| < convThreshold.
         * @param newCoords The new coordinates computed in the current round of the algorithm.
         * @param oldCoords The coordinates from the previous round of the algorithm.
         * @return @code True when converged, otherwise \endcode false.
         */
        bool isConverged(const CoordinateVector& newCoords, const CoordinateVector& oldCoords);

        /**
         * Create a weighted Laplacian matrix from G and setup the solver for this matrix.
         */
        void setupWeightedLaplacianMatrix();

        /**
         * Computes the vector L_{w,d}*x where x is stored in @a coordinates and stores the result in @a rhs.
         * @param coordinates The coordinate vector (x in the thesis)
         * @param rhs The right-hand side that stores the result of the matrix-vector multiplication.
         */
        void computeCoordinateLaplacianTerm(const CoordinateVector& coordinates, CoordinateVector& rhs);

        /**
         * Computes the repulsive forces according to Equation (8) in Gansner et al.
         * @param coordinates The current coordinates of the vertices.
         * @param b Repulsive force vector to compute
         */
        CoordinateVector computeRepulsiveForces(const CoordinateVector& coordinates, CoordinateVector& b) const;

        /**
         * Approximates the repulsive forcse by means of an octree (Barnes and Hut).
         * @param coordinates The current coordinates of the vertices.
         * @param octree Octree for Barnes-Hut approximation
         * @param theta Parameter for Barnes-Hut cell-opening criterion
         * @param b Repulsive force vector to compute
         */
        void approxRepulsiveForces(const CoordinateVector &coordinates,
                                   const Octree<double> &octree, double theta,
                                   CoordinateVector &b) const;

        /**
         * Initializes the @a coordinates corresponding to vertices in the Graph to a random point in d-dimensional space 2000^d pixel.
         * @param coordinates The coordinates of the graph.
         */
        void randomInitCoordinates(CoordinateVector& coordinates) const;

        /**
         * Initialized the @a coordinates corresponding to vertices in the Graph acoording to the random sphere placement algorithm.
         * @param coordinates The coordinates of the graph vertices.
         */
        void randomSphereCoordinates(CoordinateVector& coordinates) const;

        /**
         * Computes the squared distance ||c_i - c_j||^2 between @a coordinates @a i and @a j
         * @param coordinates
         * @param i
         * @param j
         */
        double squaredDistance(const CoordinateVector &coordinates, index i, index j) const;

        /**
         * Computes the distance ||c_i - c_j|| between @a coordinates @a i and @a j.
         * @param coordinates
         * @param i
         * @param j
         */
        inline double distance(const CoordinateVector& coordinates, const index i, const index j) const {
            return sqrt(squaredDistance(coordinates, i, j));
        }

        /**
         * Computes the squared distance ||c1_i - c2_j||^2 between coordinate @a i from @a coordinates1 and coordinate @a j from @a coordinates 2.
         * @param coordinates1
         * @param coordinates2
         * @param i
         * @param j
         */
        double squaredDistance(const CoordinateVector &coordinates1,
                               const CoordinateVector &coordinates2, index i, index j) const;

        /**
         * Computes the squared distance ||c1_i - c2_j|| between coordinate @a i from @a coordinates1 and coordinate @a j from @a coordinates 2.
         * @param coordinates1
         * @param coordinates2
         * @param i
         * @param j
         */
        inline double distance(const CoordinateVector& coordinates1, const CoordinateVector& coordinates2, const index i, const index j) const {
            return sqrt(squaredDistance(coordinates1, coordinates2, i, j));
        }

        /**
         * Computes the squared length ||c_i||^2 of coordinate @a i in @a coordinates.
         * @param coordinates
         * @param i
         */
        double squaredLength(const CoordinateVector &coordinates, index i) const;

        /**
         * Computes the length ||c_i|| of coordinate @a i in @a coordinates.
         * @param coordinates
         * @param i
         */
        inline double length(const CoordinateVector& coordinates, const index i) const {
            return sqrt(squaredLength(coordinates, i));
        }



        /**
         * Weighting factor, Gansner et al. propose 1/(edgeWeight^2).
         * @param edgeWeight The graph-theoretic distance between two vertices (i.e. the edge weight).
         * @return The corresponding weighting factor.
         */
        inline double weightingFactor(double edgeWeight) const {
            return 1.0/(edgeWeight * edgeWeight);
        }

        /**
         * Returns the sign of @a value. The sign of 0.0 is defined to be positive.
         * @param value
         */
        inline double sign(const double value) const {
            return (value >= 0.0) - (value < 0.0);
        }

        /**
         * Returns the point at index @a i stored in @a coordinates.
         * @param coordinates
         * @param i
         */
        inline Point<double> getPoint(const CoordinateVector& coordinates, index i) const {
            Point<double> p(coordinates.size());
            for (index d = 0; d < p.getDimensions(); ++d) {
                assert(coordinates[d].getDimension() > i);
                p[d] = coordinates[d][i];
            }

            return p;
        }

        /**
         * Computes the distances from each vertex to its k-hop neighborhood according to the @a graphDistance
         * type (edge weights or algebraic distance).
         * @param k
         * @param graphDistance
         */
        void computeKnownDistances(count k, GraphDistance graphDistance);

        /**
         * Add the k-neighborhood (except for the 1-neighbors) of vertex @a u to the knownDistances of @a u
         * @param u
         */
        void addKNeighborhoodOfVertex(node u, count k);

        /**
         * Computes the algebraic distances from each vertex to its k-hop neighborhood.
         * @param graph
         * @param k
         */
        void computeAlgebraicDistances(const Graph &graph, count k);
    };

} /* namespace NetworKit */
#endif // NETWORKIT_VIZ_MAXENT_STRESS_HPP_
