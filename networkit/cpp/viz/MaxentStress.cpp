// no-networkit-format
/*
 * MaxentStress.cpp
 *
 *  Created on: 22.01.2014
 *      Author: Henning Meyerhenke and Michael Wegner
 */

#include <cmath>
#include <limits>
#include <queue>

#include <networkit/auxiliary/Log.hpp>
#include <networkit/auxiliary/PrioQueue.hpp>
#include <networkit/auxiliary/Timer.hpp>
#include <networkit/components/ConnectedComponents.hpp>
#include <networkit/distance/BFS.hpp>
#include <networkit/distance/Dijkstra.hpp>
#include <networkit/numerics/LAMG/Lamg.hpp>
#include <networkit/numerics/ConjugateGradient.hpp>
#include <networkit/numerics/Preconditioner/DiagonalPreconditioner.hpp>
#include <networkit/numerics/Preconditioner/IdentityPreconditioner.hpp>
#include <networkit/viz/MaxentStress.hpp>

namespace NetworKit {

    Lamg<CSRMatrix>noneGiven(0.001);

MaxentStress::MaxentStress(const Graph &G, const count dim, const count k, double tolerance, LinearSolverType linearSolverType, bool fastComputation, GraphDistance graphDistance) : GraphLayoutAlgorithm<double>(G, dim), solver(noneGiven), resultStats(), q(0.0), alpha(1.0), alphaReduction(0.3), finalAlpha(0.008), convThreshold(0.001*0.001), coordinatesProvided(false), fastComputation(fastComputation), maxSolvesPerAlpha(50), knownDistances(std::vector<std::vector<ForwardEdge>>(G.numberOfNodes())), knownDistancesCardinality(0), dim(dim), hasRun(false) {

    switch (linearSolverType) {
        case LAMG:
            this->solver = Lamg<CSRMatrix>(tolerance);
            break;
        case CONJUGATE_GRADIENT_IDENTITY_PRECONDITIONER:
            this->solver = ConjugateGradient<CSRMatrix, IdentityPreconditioner>(tolerance);
            break;
        case CONJUGATE_GRADIENT_DIAGONAL_PRECONDITIONER:
            this->solver = ConjugateGradient<CSRMatrix, DiagonalPreconditioner>(tolerance);
    }
    computeKnownDistances(k, graphDistance);
}

MaxentStress::MaxentStress(const Graph &G, const count dim, const std::vector<Point<double>> &coordinates, const count k, double tolerance, LinearSolverType linearSolverType, bool fastComputation, GraphDistance graphDistance) : GraphLayoutAlgorithm<double>(G, dim), solver(noneGiven), resultStats(), q(0.0), alpha(1.0), alphaReduction(0.3), finalAlpha(0.008), convThreshold(0.001*0.001), coordinatesProvided(true), fastComputation(fastComputation), maxSolvesPerAlpha(50), knownDistances(std::vector<std::vector<ForwardEdge>>(G.numberOfNodes())), knownDistancesCardinality(0), dim(dim), hasRun(false) {
    switch (linearSolverType) {
        case LAMG:
            this->solver = Lamg<CSRMatrix>(tolerance);
            break;
        case CONJUGATE_GRADIENT_IDENTITY_PRECONDITIONER:
            this->solver = ConjugateGradient<CSRMatrix, IdentityPreconditioner>(tolerance);
            break;
        case CONJUGATE_GRADIENT_DIAGONAL_PRECONDITIONER:
            this->solver = ConjugateGradient<CSRMatrix, DiagonalPreconditioner>(tolerance);
    }
    vertexCoordinates = coordinates;
    computeKnownDistances(k, graphDistance);
}

MaxentStress::ResultStats MaxentStress::runAlgo() {
    run();
    return this->resultStats;
}

void MaxentStress::run() {
    // Check if the graph is connected. We currently can't handle unconnected graphs.
    ConnectedComponents cc(*G);
    cc.run();
    if (cc.numberOfComponents() != 1) {
        throw std::invalid_argument( "ERROR: The supplied graph is not connected. Currently MaxentStress only handles connected graphs.");
    }

    Aux::Timer t;
    double solveTime = 0;
    double rhsTime = 0;
    double approxTime = 0.0;
    t.start();
    setupWeightedLaplacianMatrix(); // create weighted Laplacian matrix and setup the solver
    t.stop();
    solveTime += t.elapsedMicroseconds();
    CoordinateVector oldCoordinates(dim, Vector(G->upperNodeIdBound()));

    if (!coordinatesProvided && !hasRun) {
        randomSphereCoordinates(oldCoordinates);
    } else {
#pragma omp parallel for
        for (omp_index i = 0; i < static_cast<omp_index>(vertexCoordinates.size()); ++i) {
            for (index d = 0; d < dim; ++d) {
                oldCoordinates[d][i] = vertexCoordinates[i][d];
            }
        }
    }

    CoordinateVector newCoordinates = oldCoordinates;

    Aux::Timer timer;
    timer.start();
    double currentAlpha = alpha;
    bool converged = false;

    CoordinateVector repulsiveForces(dim, Vector(G->numberOfNodes(), 0));
    count currentLowerBound = 0;
    count newLowerBound = 0;
    while (!converged) { // solve up to maxSolvesPerAlpha linear systems
        INFO("Running with alpha = ", currentAlpha);
        for (count numSolves = 0; numSolves < maxSolvesPerAlpha; ++numSolves) {
            oldCoordinates = newCoordinates;

            t.start();
            newLowerBound = floor(5 * std::log(numSolves));
            if (newLowerBound != currentLowerBound) {
                repulsiveForces = CoordinateVector(dim, Vector(G->numberOfNodes(), 0));
                Octree<double> octree(oldCoordinates);
                approxRepulsiveForces(oldCoordinates, octree, 0.6, repulsiveForces);
                currentLowerBound = newLowerBound;
            }
            t.stop();
            approxTime += t.elapsedMicroseconds();

            t.start();
            CoordinateVector rhs(dim, Vector(G->numberOfNodes()));
            computeCoordinateLaplacianTerm(oldCoordinates, rhs);
            t.stop();
            rhsTime += t.elapsedMicroseconds();

            t.start();
            for (index d = 0; d < dim; ++d) {
                if (numSolves < maxSolvesPerAlpha/5) {
                    rhs[d] /= rhs[d].length();
                }
                rhs[d] += currentAlpha * repulsiveForces[d];
            }
            t.stop();
            rhsTime += t.elapsedMicroseconds();

            // correcting rhs to be zero-sum
            Point<double> sum(oldCoordinates.size());
            for (index d = 0; d < dim; ++d) {
                for (index i = 0; i < G->numberOfNodes(); ++i) {
                    sum[d] += rhs[d][i];
                }
                sum[d] /= G->numberOfNodes();
            }

#pragma omp parallel for
            for (omp_index i = 0; i < static_cast<omp_index>(G->numberOfNodes()); ++i) {
                for (index d = 0; d < dim; ++d) {
                    rhs[d][i] -= sum[d];
                }
            }

            t.start();
            solver.parallelSolve(rhs, newCoordinates, 300000, 30);
            t.stop();
            solveTime += t.elapsedMicroseconds();
            converged = isConverged(newCoordinates, oldCoordinates);
            if (converged) {
                if (!fastComputation) {
                    converged = false;
                } else {
                    INFO("Converged after ", numSolves+1, " solves");
                    break;
                }
            }
        }

        INFO("converged: ", converged);

        currentAlpha *= alphaReduction; // Cooling: reduce alpha for next round
        converged = converged || currentAlpha < finalAlpha;
    }

    // write coordinates to vertexCoordinates vector
#pragma omp parallel for
    for (omp_index i = 0; i < static_cast<omp_index>(G->upperNodeIdBound()); ++i) {
        for (index d = 0; d < dim; ++d) {
            this->vertexCoordinates[i][d] = newCoordinates[d][i];
        }
    }
    timer.stop();

    hasRun = true;

    this->resultStats.rhsTime = rhsTime / 1000;
    this->resultStats.solveTime = solveTime / 1000;
    this->resultStats.approxEntropyTerm = approxTime / 1000;

    INFO("Time spent on rhs ", rhsTime/1000, " and solve time was ", solveTime/1000);
    INFO("approxTime: ", approxTime/1000);
    INFO("Graph drawn in ", timer.elapsedMilliseconds());
}

double MaxentStress::computeScalingFactor() {
    count n = G->numberOfNodes();
    bool weighted = false;
    Graph augmentedGraph(n, true, false);
    for (node u = 0; u < n; ++u) {
        for (ForwardEdge &edge : knownDistances[u]) {
            augmentedGraph.addEdge(u, edge.head, edge.weight);
            weighted = weighted || edge.weight != 1.0;
        }
    }

    double topFraction = 0.0;
    double bottomFraction = 0.0;
#pragma omp parallel for reduction(+:topFraction)
    for (omp_index u = 0; u < static_cast<omp_index>(n); ++u) {
        std::unique_ptr<SSSP> sssp = weighted? std::move(std::unique_ptr<SSSP>(new Dijkstra(augmentedGraph, u, false, false))) : std::move(std::unique_ptr<SSSP>(new BFS(augmentedGraph, u, false, false)));
        sssp->run();
        augmentedGraph.forNodes([&](node v) {
            double geometricDist = this->vertexCoordinates[u].distance(this->vertexCoordinates[v]);
            if (sssp->distance(v) < 1e-5) return;
            topFraction += geometricDist / sssp->distance(v);
        });
    }

#pragma omp parallel for reduction(+:bottomFraction)
    for (omp_index u = 0; u < static_cast<omp_index>(n); ++u) {
        std::unique_ptr<SSSP> sssp = weighted? std::move(std::unique_ptr<SSSP>(new Dijkstra(augmentedGraph, u, false, false))) : std::move(std::unique_ptr<SSSP>(new BFS(augmentedGraph, u, false, false)));
        sssp->run();
        augmentedGraph.forNodes([&](node v) {
            double squaredGeometricDist = this->vertexCoordinates[u].squaredDistance(this->vertexCoordinates[v]);
            if (sssp->distance(v) < 1e-5) return;
            bottomFraction += squaredGeometricDist / (sssp->distance(v) * sssp->distance(v));
        });
    }
    return topFraction / bottomFraction;
}

void MaxentStress::scaleLayout() {
    double s = computeScalingFactor();
    INFO("Scaling factor: ", s);

    // scale vertex coordinates with s
#pragma omp parallel for
    for (omp_index i = 0; i < static_cast<omp_index>(this->vertexCoordinates.size()); ++i) {
        this->vertexCoordinates[i].scale(s);
    }
}

double MaxentStress::fullStressMeasure() {
    double energy = 0.0;

    count n = G->numberOfNodes();
    bool weighted = false;
    Graph augmentedGraph(n, true, false);
    for (node u = 0; u < n; ++u) {
        for (ForwardEdge &edge : knownDistances[u]) {
            augmentedGraph.addEdge(u, edge.head, edge.weight);
            weighted = weighted || edge.weight != 1.0;
        }
    }

#pragma omp parallel for reduction(+:energy)
    for (omp_index u = 0; u < static_cast<omp_index>(n); ++u) {
        std::unique_ptr<SSSP> sssp = weighted? std::move(std::unique_ptr<SSSP>(new Dijkstra(augmentedGraph, u, false, false))) : std::move(std::unique_ptr<SSSP>(new BFS(augmentedGraph, u, false, false)));
        sssp->run();
        G->forNodes([&](node v) {
            double geometricDist = this->vertexCoordinates[u].distance(this->vertexCoordinates[v]);
            if (sssp->distance(v) < 1e-5) return;
            energy += (geometricDist - sssp->distance(v))*(geometricDist - sssp->distance(v))/(sssp->distance(v)*sssp->distance(v));
        });
    }

    return energy/2;
}

double MaxentStress::maxentMeasure() {
    double energy  = 0;
    double entropy = 0;
    count n = G->numberOfNodes();

    Graph augmentedGraph(n, true, false);
    for (node u = 0; u < n; ++u) {
        for (ForwardEdge &edge : knownDistances[u]) {
            augmentedGraph.addEdge(u, edge.head, edge.weight);
        }
    }

#pragma omp parallel for reduction(+:entropy)
    for (omp_index u = 0; u < static_cast<omp_index>(n); ++u) {
        augmentedGraph.forNodes([&](node v) {
            if (u == v) return;
            double dist = std::max(this->vertexCoordinates[u].distance(this->vertexCoordinates[v]), 1e-5);
            entropy += (fabs(q) < 0.001)? log(dist) : pow(dist, -q);
        });
    }

    INFO("entropy: ", entropy);

    for (node u = 0; u < n; ++u) {
        augmentedGraph.forNeighborsOf(u, [&](node v, edgeweight w) {
            double dist = std::max(this->vertexCoordinates[u].distance(this->vertexCoordinates[v]), 1e-5);
            energy += (dist - w)*(dist - w)/(w*w);
            entropy -= (fabs(q) < 0.001)? log(dist) : pow(dist, -q);
        });
    }

    if(fabs(q) > 0.001) {
        entropy *= -sign(q);
    }

    energy -= finalAlpha*entropy;
    return energy/2;
}

double MaxentStress::meanDistanceError() {
    double sum = 0.0;
    for (node u = 0; u < knownDistances.size(); ++u) {
        for (const ForwardEdge& edge : knownDistances[u]) {
            sum += fabs(vertexCoordinates[u].distance(vertexCoordinates[edge.head]) - edge.weight) / edge.weight;
        }
    }

    return sum / knownDistancesCardinality;
}

double MaxentStress::ldme() {
    double sum = 0.0;
    for (node u = 0; u < knownDistances.size(); ++u) {
        for (const ForwardEdge& edge : knownDistances[u]) {
            double coordDist = vertexCoordinates[u].distance(vertexCoordinates[edge.head]);
            sum += (edge.weight - coordDist) * (edge.weight - coordDist);
        }
    }

    sum /= knownDistancesCardinality;
    return sqrt(sum);
}

bool MaxentStress::isConverged(const CoordinateVector& newCoords, const CoordinateVector& oldCoords) {
    assert(newCoords.size() == oldCoords.size());

    double relChange = 0.0;
    double oldCoordsSqLength = 0.0;
#pragma omp parallel for reduction(+:relChange,oldCoordsSqLength)
    for (omp_index i = 0; i < static_cast<omp_index>(newCoords[0].getDimension()); ++i) {
        relChange += squaredDistance(newCoords, oldCoords, i, i);
        oldCoordsSqLength += squaredLength(oldCoords, i);
    }

    return relChange / oldCoordsSqLength < convThreshold;
}


void MaxentStress::setupWeightedLaplacianMatrix() {
    count n = G->numberOfNodes();
    std::vector<index> rowIdx(n+1, 0);
    std::vector<index> columnIdx(n + knownDistancesCardinality, 0); // currently only supports simple graphs
    std::vector<double> nonzeros(columnIdx.size());

    index idx = 0;
    for (index i = 0; i < n; ++i) {
        double weightedDegree = 0.0;
        for (ForwardEdge &edge : knownDistances[i]) {
            double weightFactor = weightingFactor(edge.weight);

            columnIdx[idx] = edge.head;
            nonzeros[idx] = -weightFactor;
            weightedDegree += weightFactor;
            idx++;
        }

        // add diagonal element
        columnIdx[idx] = i;
        nonzeros[idx] = weightedDegree;
        idx++;

        rowIdx[i+1] = knownDistances[i].size()+1; // +1 for diagonal
    }

    // compute correct rowIdx offsets
    for (index i = 1; i < rowIdx.size(); ++i) {
        rowIdx[i] += rowIdx[i-1];
    }

    CSRMatrix laplacian(n, n, rowIdx, columnIdx, nonzeros);
    solver.setupConnected(laplacian);
}

void MaxentStress::computeCoordinateLaplacianTerm(const CoordinateVector& coordinates, CoordinateVector& rhs) {
    count n = G->numberOfNodes();
#pragma omp parallel for
    for (omp_index i = 0; i < static_cast<omp_index>(n); ++i) {
        double weightedDegree = 0.0;
        for (const ForwardEdge& edge : knownDistances[i]) {
            double dist = std::max(distance(coordinates, i, edge.head), 1e-5);
            double w = weightingFactor(edge.weight) * edge.weight / dist; // w_{ij} * d_{i,j} / ||x_i - x_j|| NOTE: The last term is multiplied in the Gansner et al. paper which is wrong!
            for (index d = 0; d < dim; ++d) {
                rhs[d][i] += -w * coordinates[d][edge.head];
            }
            weightedDegree += w;
        }

        for (index d = 0; d < dim; ++d) {
            rhs[d][i] += weightedDegree * coordinates[d][i];
        }
    }
}

CoordinateVector MaxentStress::computeRepulsiveForces(const CoordinateVector& coordinates, CoordinateVector &b) const {
    count n = G->numberOfNodes();
    double qSign = sign(this->q);
    double q2 = (q+2)/2;

#pragma omp parallel for
    for (omp_index i = 0; i < static_cast<omp_index>(n); ++i) {
        std::vector<bool> knownDist(n, false);
        for (const ForwardEdge &edge : knownDistances[i]) {
            knownDist[edge.head] = true;
        }

        for (index j = 0; j < n; ++j) {
            if (!knownDist[j] && i != j) {
                double sqDist = std::max(squaredDistance(coordinates, i, j), 1e-3); // ||x_i - x_j||
                double factor = qSign * 1.0/std::pow(sqDist, q2);
                for (index d = 0; d < dim; ++d) {
                    b[d][i] += factor * (coordinates[d][i] - coordinates[d][j]); // sum_{\{i,k\} \in S} dist^{-q-2} * (x_{i,d} - x_{j,d})
                }
            }
        }
    }

    // normalize b
    for (index d = 0; d < dim; ++d) {
        b[d] /= b[d].length();
    }

    return b;
}

void MaxentStress::approxRepulsiveForces(const CoordinateVector& coordinates, const Octree<double>& octree, const double theta, CoordinateVector& b) const {
    count n = G->numberOfNodes();
    double qSign = sign(q);
    double q2 = (q+2)/2;

#pragma omp parallel for
    for (omp_index i = 0; i < static_cast<omp_index>(n); ++i) {
        Point<double> pI = getPoint(coordinates, i);
        auto approximateNeighbor = [&](const count numNodes, const Point<double>& centerOfMass, const double sqDist) {
            if (sqDist < 1e-5) return;
            double factor = qSign * numNodes * 1.0/pow(sqDist, q2);
            for (index d = 0; d < dim; ++d) {
                b[d][i] += factor * (pI[d] - centerOfMass[d]);
            }
        };

        octree.approximateDistance(pI, theta, approximateNeighbor);
    }

    // normalize b
    for (index d = 0; d < dim; ++d) {
        b[d] /= b[d].length();
    }
}

void MaxentStress::computeKnownDistances(const count k, const GraphDistance graphDistance) {
    switch (graphDistance) {
    case EDGE_WEIGHT:
        // add edges to known distances
        G->parallelForNodes([&](node u) {
            G->forNeighborsOf(u, [&](node v, edgeweight w) {
                knownDistances[u].push_back({v, w});
            });
            // add k-neighborhood
            if (k > 1) { // 1-neighborhood are adjacent vertices that are already added above
                addKNeighborhoodOfVertex(u, k);
            }
        });
        break;
    case ALGEBRAIC_DISTANCE:
        computeAlgebraicDistances(*G, k);
        break;
    default:
        break;
    }

    // finally determine cardinality of knownDistances
    count cardinality = 0;
#pragma omp parallel for reduction(+:cardinality)
    for (omp_index i = 0; i < static_cast<omp_index>(knownDistances.size()); ++i) {
        cardinality += knownDistances[i].size();
    }

    knownDistancesCardinality = cardinality;
    INFO("|S| = ", cardinality);

    // check if graph has more than 30 percent degree-1 vertices: if yes, set q to 0.8
    count numDegreeOneNodes = 0;
    G->forNodes([&](node u) {
        numDegreeOneNodes += G->degree(u) == 1;
    });

    if ((double) numDegreeOneNodes / (double) G->numberOfNodes() > 0.3) {
        q = 0.8;
        INFO("Setting q to 0.8 because we have ", (double) numDegreeOneNodes / (double) G->numberOfNodes() * 100, " % degree-1 nodes");
    }
}

void MaxentStress::addKNeighborhoodOfVertex(const node u, const count k) {
    if (G->isWeighted()) {
        // first compute depths with bfs
        std::queue<node> Q;
        std::vector<count> depth(G->numberOfNodes(), none);
        count numVerticesToVisit = 0;

        Q.push(u);
        depth[u] = 0;
        while (!Q.empty()) {
            node v = Q.front(); Q.pop();
            G->forNeighborsOf(v, [&](node w, edgeweight) {
                if (depth[w] == none) { // node has not been visited yet
                    depth[w] = depth[v]+1;
                    if (depth[w] > 1) numVerticesToVisit++;

                    if (depth[w] < k) { // only push neighbors that are less than k steps away
                        Q.push(w);
                    }
                }
            });
        }

        Aux::PrioQueue<edgeweight, node> PQ(G->numberOfNodes());
        std::vector<edgeweight> dist(G->numberOfNodes(), std::numeric_limits<edgeweight>::max());

        dist[u] = 0;
        PQ.insert(0, u);
        while (!PQ.empty() && numVerticesToVisit > 0) {
            node v = PQ.extractMin().second;
            if (1 < depth[v] && depth[v] <= k) {
                knownDistances[u].push_back({v, dist[v]});
                numVerticesToVisit--;
            }
            G->forNeighborsOf(v, [&](node w, edgeweight weight) {
                if (dist[v] + weight < dist[w]) {
                    dist[w] = dist[v] + weight;
                    PQ.changeKey(dist[w], w);
                }
            });
        }
    } else {
        std::queue<node> Q;
        std::vector<count> depth(this->G->numberOfNodes(), none);

        Q.push(u);
        depth[u] = 0;
        while (!Q.empty()) {
            node v = Q.front(); Q.pop();
            G->forNeighborsOf(v, [&](node w, edgeweight) {
                if (depth[w] == none) { // node has not been visited yet
                    depth[w] = depth[v]+1;

                    if (depth[w] < k) { // only push neighbors that are less than k steps away
                        Q.push(w);
                    }

                    if (1 < depth[w] && depth[w] <= k) { // add neighbors that are more than 1 and less than or equal to k steps away
                        knownDistances[u].push_back({w, (double)depth[w]});
                    }
                }
            });
        }
    }
}

void MaxentStress::computeAlgebraicDistances(const Graph& graph, const count k) {
    AlgebraicDistance ad(graph);
    ad.preprocess();

    std::vector<double> minDist(G->numberOfNodes(), 1.0);
    std::vector<double> maxDist(G->numberOfNodes(), 0.0);
    std::vector<double> distances;
    graph.parallelForNodes([&](node u) {
        std::queue<node> Q;
        std::vector<count> depth(G->numberOfNodes(), none);

        Q.push(u);
        depth[u] = 0;
        while (!Q.empty()) {
            node v = Q.front(); Q.pop();
            graph.forNeighborsOf(v, [&](node w, edgeweight) {
                if (depth[w] == none) { // node has not been visited yet
                    depth[w] = depth[v]+1;

                    if (depth[w] < k) { // only push neighbors that are less than k steps away
                        Q.push(w);
                    }

                    if (depth[w] <= k) { // add neighbors that are less than or equal to k steps away
                        double algebraicDist = ad.distance(u,w);
                        if (algebraicDist == 0.0) {
                            algebraicDist = 1e-5;
                        }
                        algebraicDist /= sqrt(static_cast<double>(G->degree(u) * G->degree(w)));
                        knownDistances[u].push_back({w, algebraicDist});
                        if (std::isnan(algebraicDist)) INFO("Warning: nan dist");
                        minDist[u] = std::min(minDist[u], algebraicDist);
                        maxDist[u] = std::max(maxDist[u], algebraicDist);
                    }
                }
            });
        }
    });

    double minimumDist = minDist[0];
    double maximumDist = maxDist[0];
    for (node u = 1; u < G->numberOfNodes(); ++u) {
        minimumDist = std::min(minimumDist, minDist[u]);
        maximumDist = std::max(maximumDist, maxDist[u]);
    }

    INFO("[min, max] = [", minimumDist, ",", maximumDist, "]");

    graph.parallelForNodes([&](node u) {
        for (index i = 0; i < knownDistances[u].size(); ++i) {
            knownDistances[u][i].weight = std::log(2.0 + (knownDistances[u][i].weight - minimumDist) / (maximumDist-minimumDist) * 11);
        }
    });

    for (node u = 0; u < G->numberOfNodes(); ++u) {
        for (const ForwardEdge& edge : knownDistances[u]) {
            node v = edge.head;
            if (u < v) {
                bool backEdgePresent = false;
                for (const ForwardEdge& bEdge : knownDistances[v]) {
                    if (bEdge.head == u) {
                        backEdgePresent = true;
                        break;
                    }
                }
                if (!backEdgePresent) INFO("WARNING: Missing backEdge for edge (", u, ",", v, ")");
            }

        }
    }
}

void MaxentStress::randomInitCoordinates(CoordinateVector &coordinates) const {

#pragma omp parallel for
    for (omp_index i = 0; i < static_cast<omp_index>(coordinates[0].getDimension()); ++i) {
        for (index d = 0; d < dim; ++d) {
            coordinates[d][i] = Aux::Random::real() * 50; // 50 x 50 pixel
        }
    }
}

void MaxentStress::randomSphereCoordinates(CoordinateVector &coordinates) const {
    // find node with highest degree
    node maxDegNode = 0;
    count maxDeg = G->degree(maxDegNode);
    G->forNodes([&](node u) {
        if (G->degree(u) > maxDeg) {
            maxDegNode = u;
            maxDeg = G->degree(u);
        }
    });


    // set coordinate of node 0 to (0,0,...,0)
    for (index d = 0; d < dim; ++d) {
        coordinates[d][maxDegNode] = 0.0;
    }

    std::vector<bool> coordinateSet(G->upperNodeIdBound(), false);
    coordinateSet[maxDegNode] = true;
    count numSet = 1;

    while (numSet < G->numberOfNodes()) {
        node start = none;
        G->forNodes([&](node u) {
            if (coordinateSet[u]) {
                start = u;
                return;
            }
        });

        // perform BFS from start
        std::queue<node> Q;
        Q.push(start);
        while (!Q.empty()) {
            node u = Q.front(); Q.pop();
            G->forNeighborsOf(u, [&](node u, node v, edgeweight w, index) {
                if (!coordinateSet[v]) {
                    Vector p(dim);
                    for (index d = 0; d < dim; ++d) {
                        p[d] = 2*Aux::Random::real() - 1;
                    }
                    p *= w / p.length();
                    for (index d = 0; d < dim; ++d) {
                        coordinates[d][v] = coordinates[d][u] + p[d];
                    }
                    coordinateSet[v] = true;
                    numSet++;
                    Q.push(v);
                }
            });
        }
    }
}

double MaxentStress::squaredDistance(const CoordinateVector& coordinates, const index i, const index j) const {
    double dist = 0.0;
    for (index d = 0; d < dim; ++d) {
        double diff = coordinates[d][i] - coordinates[d][j];
        dist += diff * diff;
    }
    return dist;
}

double MaxentStress::squaredDistance(const CoordinateVector& coordinates1, const CoordinateVector& coordinates2, const index i, const index j) const {
    double dist = 0.0;
    for (index d = 0; d < dim; ++d) {
        double diff = coordinates1[d][i] - coordinates2[d][j];
        dist += diff * diff;
    }

    return dist;
}

double MaxentStress::squaredLength(const CoordinateVector& coordinates, const index i) const {
    double length = 0.0;
    for (index d = 0; d < dim; ++d) {
        length += coordinates[d][i] * coordinates[d][i];
    }

    return length;
}
} /* namespace NetworKit */
