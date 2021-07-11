// no-networkit-format
/*
 * AlgebraicDistance.cpp
 *
 *  Created on: 03.11.2015
 *      Author: Henning Meyerhenke, Christian Staudt, Michael Hamann
 */

#include <omp.h>

#include <networkit/auxiliary/Timer.hpp>
#include <networkit/distance/AlgebraicDistance.hpp>


namespace NetworKit {

AlgebraicDistance::AlgebraicDistance(const Graph& G, count numberSystems, count numberIterations, double omega, index norm, bool withEdgeScores) : NodeDistance(G), numSystems(numberSystems), numIters(numberIterations), omega(omega), norm(norm), withEdgeScores(withEdgeScores) {
    if ((omega < 0.0) || (omega > 1.0)) throw std::invalid_argument("omega must be in [0,1]");
    if (withEdgeScores && !G.hasEdgeIds()) {
        throw std::runtime_error("edges have not been indexed - call indexEdges first");
    }
}

void AlgebraicDistance::randomInit() {
    // allocate space for loads
    loads.resize(numSystems*G->upperNodeIdBound());

    #pragma omp parallel for
    for (omp_index i = 0; i < static_cast<omp_index>(loads.size()); ++i) {
        loads[i] = Aux::Random::real();
    }
}

void AlgebraicDistance::preprocess() {
    Aux::Timer running1;
    running1.start();
    // random init
    randomInit();

    // main loop

    {
        std::vector<double> oldLoads(loads.size());

        for (index iter = 0; iter < numIters; ++iter) {
            // store previous iteration
            loads.swap(oldLoads);

            G->balancedParallelForNodes([&](node u) {
                std::vector<double> val(numSystems, 0.0);

                double weightedDeg = 0;
                // step 1
                G->forNeighborsOf(u, [&](node v, edgeweight weight) {
                    for (index i = 0; i < numSystems; ++i) {
                        val[i] += weight * oldLoads[v*numSystems + i];
                    }

                    weightedDeg += weight;
                });

                for (index i = 0; i < numSystems; ++i) {
                    val[i] /= weightedDeg;

                    // step 2
                    loads[u*numSystems + i] = (1 - omega) * oldLoads[u*numSystems + i] + omega * val[i];
                }
            });
        }
    }

    // normalization. Compute min/max over all nodes per system (and per thread)
    std::vector<std::vector<double>> minPerThread(omp_get_max_threads(), std::vector<double>(numSystems, std::numeric_limits<double>::max()));
    std::vector<std::vector<double>> maxPerThread(omp_get_max_threads(), std::vector<double>(numSystems, std::numeric_limits<double>::lowest()));
    G->parallelForNodes([&](node u) {
        auto tid = omp_get_thread_num();
        const index startId = u*numSystems;
        for (index sys = 0; sys < numSystems; ++sys) {
            minPerThread[tid][sys] = std::min(minPerThread[tid][sys], loads[startId + sys]);
            maxPerThread[tid][sys] = std::max(maxPerThread[tid][sys], loads[startId + sys]);
        }
    });

    std::vector<double> minPerSystem = std::move(minPerThread[0]);
    std::vector<double> maxPerSystem = std::move(maxPerThread[0]);
    for (index i = 1; i < minPerThread.size(); ++i) {
        for (index sys = 0; sys < numSystems; ++sys) {
            minPerSystem[sys] = std::min(minPerSystem[sys], minPerThread[i][sys]);
            maxPerSystem[sys] = std::max(maxPerSystem[sys], maxPerThread[i][sys]);
        }
    }

    // set normalized values: new = (min - old) / (min - max)
    // normalization is per system
    G->parallelForNodes([&](node u) {
        const index startId = u*numSystems;
        for (index sys = 0; sys < numSystems; ++sys) {
            loads[startId + sys] = (minPerSystem[sys] - loads[startId + sys]) / (minPerSystem[sys] - maxPerSystem[sys]);
        }
    });

    // calculate edge scores

    if (withEdgeScores) {
        edgeScores.resize(G->upperEdgeIdBound(), std::numeric_limits<edgeweight>::max());

        G->parallelForEdges([&](node u, node v, edgeid eid) {
            edgeScores[eid] = distance(u, v);
        });
    }


    running1.stop();
    INFO("elapsed millisecs for AD preprocessing: ", running1.elapsedMilliseconds(), "\n");
}

double AlgebraicDistance::distance(node u, node v) {
    if (loads.empty()) {
        throw std::runtime_error("Call preprocess() first.");
    }
    double result = 0.0;

    if (norm == MAX_NORM) {
        for (index sys = 0; sys < numSystems; ++sys) {
            double absDiff = fabs(loads[u*numSystems + sys] - loads[v*numSystems + sys]);
            if (absDiff > result) {
                result = absDiff;
            }
        }
    } else {
        for (index sys = 0; sys < numSystems; ++sys) {
            double absDiff = fabs(loads[u*numSystems + sys] - loads[v*numSystems + sys]);
            result += pow(absDiff, norm);
        }
        result = pow(result, 1.0 / (double) norm);
    }

    return std::isnan(result) ? 0 : result;
}

std::vector<double> AlgebraicDistance::getEdgeScores() {
    if (!withEdgeScores) throw std::runtime_error("set constructor parameter 'withEdgeScores' to true");
    return edgeScores;
}

} /* namespace NetworKit */
