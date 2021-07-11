// no-networkit-format
/*
 * HyperbolicGenerator.cpp
 *
 *      Authors: Mustafa zdayi and Moritz v. Looz
 *
 * This generator contains algorithms described in two publications.
 *
 * For T=0, the relevant publication is
 * "Generating massive complex networks with hyperbolic geometry faster in practice" by
 * Moritz von Looz, Mustafa zdayi, Sren Laue and Henning Meyerhenke, presented at HPEC 2016.
 *
 * For T>0, it is
 * "Querying Probabilistic Neighborhoods in Spatial Data Sets Efficiently" by Moritz von Looz
 * and Henning Meyerhenke, presented at IWOCA 2016.
 *
 * The model of hyperbolic random graphs is presented in
 * "Hyperbolic geometry of complex networks. Physical Review E, 82:036106, Sep 2010." by
 *   Dmitri Krioukov, Fragkiskos Papadopoulos, Maksim Kitsak, Amin Vahdat, and Marian Boguna
 *
 */

#include <algorithm>
#include <omp.h>
#include <random>

#include <networkit/auxiliary/Parallel.hpp>
#include <networkit/auxiliary/Random.hpp>
#include <networkit/generators/HyperbolicGenerator.hpp>
#include <networkit/generators/quadtree/Quadtree.hpp>
#include <networkit/graph/GraphBuilder.hpp>

namespace NetworKit {

/**
 * Construct a generator for n nodes and m edges
 */
HyperbolicGenerator::HyperbolicGenerator(count n, double avgDegree, double plexp, double T) {
    nodeCount = n;
    if (plexp <= 2) throw std::runtime_error("Exponent of power-law degree distribution must be > 2");
    if (T < 0 || T == 1) throw std::runtime_error("Temperature must be non-negative and not 1.");//Really necessary? Graphs with T=1 can be generated, only their degree is not controllable
    if (avgDegree >= n) throw std::runtime_error("Average Degree must be at most n-1");
    if (T < 1) {
        alpha = 0.5*(plexp-1);
    } else {
        alpha = 0.5*(plexp-1)/T;
    }

    R = HyperbolicSpace::getTargetRadius(n, n*avgDegree/2, alpha, T);
    temperature=T;
    initialize();
}

void HyperbolicGenerator::initialize() {
    if (temperature == 0) {
        capacity = 1000;
    } else {
        capacity = 10;
    }
    theoreticalSplit = false;
    threadtimers.resize(omp_get_max_threads());
    balance = 0.5;
}

Graph HyperbolicGenerator::generate() {
    return generate(nodeCount, R, alpha, temperature);
}

Graph HyperbolicGenerator::generate(count n, double R, double alpha, double T) {
    assert(R > 0);
    vector<double> angles(n);
    vector<double> radii(n);

    //sample points randomly
    HyperbolicSpace::fillPoints(angles, radii, R, alpha);
    vector<index> permutation(n);

    index p = 0;
    std::generate(permutation.begin(), permutation.end(), [&p](){return p++;});

    //can probably be parallelized easily, but doesn't bring much benefit
    Aux::Parallel::sort(permutation.begin(), permutation.end(), [&angles,&radii](index i, index j){return angles[i] < angles[j] || (angles[i] == angles[j] && radii[i] < radii[j]);});

    vector<double> anglecopy(n);
    vector<double> radiicopy(n);

    #pragma omp parallel for
    for (omp_index j = 0; j < static_cast<omp_index>(n); j++) {
        anglecopy[j] = angles[permutation[j]];
        radiicopy[j] = radii[permutation[j]];
    }

    INFO("Generated Points");
    return generate(anglecopy, radiicopy, R, T);
}

Graph HyperbolicGenerator::generateCold(const vector<double> &angles, const vector<double> &radii, double R) {
    const count n = angles.size();
    assert(radii.size() == n);

    for (index i = 0; i < n; i++) {
        assert(radii[i] < R);
    }

    vector<index> permutation(n);
    #pragma omp parallel for
    for (omp_index i = 0; i < static_cast<omp_index>(n); i++) {
        permutation[i] = i;
    }

    //can probably be parallelized easily, but doesn't bring much benefit
    Aux::Parallel::sort(permutation.begin(), permutation.end(), [&angles,&radii](index i, index j){return angles[i] < angles[j] || (angles[i] == angles[j] && radii[i] < radii[j]);});

    vector<double> bandRadii = getBandRadii(n, R);
    //Initialize empty bands
    vector<vector<Point2DWithIndex<double>>> bands(bandRadii.size() - 1);
    //Put points to bands
    #pragma omp parallel for
    for (omp_index j = 0; j < static_cast<omp_index>(bands.size()); j++){
        for (index i = 0; i < n; i++){
            double alias = permutation[i];
            if (radii[alias] >= bandRadii[j] && radii[alias] <= bandRadii[j+1]){
                bands[j].push_back(Point2DWithIndex<double>(angles[alias], radii[alias], alias));
            }
        }
    }

    const count bandCount = bands.size();
    const double coshR = cosh(R);
    assert(radii.size() == n);

    Aux::Timer bandTimer;
    bandTimer.start();

    //1.Extract band angles to use them later, can create a band class to handle this more elegantly
    vector<vector<double>> bandAngles(bandCount);
    #pragma omp parallel for
    for (omp_index i=0; i < static_cast<omp_index>(bandCount); i++){
        const count currentBandSize = bands[i].size();
        bandAngles[i].resize(currentBandSize);
        for(index j=0; j < currentBandSize; j++) {
            bandAngles[i][j] = bands[i][j].getX();
        }
        if (!std::is_sorted(bandAngles[i].begin(), bandAngles[i].end())) {
            throw std::runtime_error("Points in bands must be sorted.");
        }
    }
    bandTimer.stop();
    INFO("Extracting band angles took ", bandTimer.elapsedMilliseconds(), " milliseconds.");

    //2.Insert edges
    Aux::Timer timer;
    timer.start();
    vector<double> empty;
    GraphBuilder result(n, false, false);

    #pragma omp parallel
    {
        index id = omp_get_thread_num();
        threadtimers[id].start();
        #pragma omp for schedule(guided) nowait
        for (omp_index i = 0; i < static_cast<omp_index>(n); i++) {
            const double coshr = cosh(radii[i]);
            const double sinhr = sinh(radii[i]);
            count expectedDegree = (4/PI)*n*exp(-(radii[i])/2);
            vector<index> near;
            near.reserve(expectedDegree*1.1);
            Point2DWithIndex<double> pointV(angles[i], radii[i], i);
            for(index j = 0; j < bandCount; j++){
                if(directSwap || bandRadii[j+1] > radii[i]){
                    double minTheta, maxTheta;
                    std::tie (minTheta, maxTheta) = getMinMaxTheta(angles[i], radii[i], bandRadii[j], R);
                    vector<Point2DWithIndex<double>> neighborCandidates = getPointsWithinAngles(minTheta, maxTheta, bands[j], bandAngles[j]);

                    const count sSize = neighborCandidates.size();
                    for(index w = 0; w < sSize; w++){
                        double deltaPhi = PI - abs(PI-abs(angles[i] - neighborCandidates[w].getX()));
                        if (coshr*cosh(neighborCandidates[w].getY())-sinhr*sinh(neighborCandidates[w].getY())*cos(deltaPhi) <= coshR) {
                            if (neighborCandidates[w].getIndex() != i){
                                near.push_back(neighborCandidates[w].getIndex());
                            }
                        }
                    }
                }
            }
            if (directSwap) {
                auto newend = std::remove(near.begin(), near.end(), i); //no self loops!
                if (newend != near.end()) {
                    assert(newend+1 == near.end());
                    assert(*(newend)==i);
                    near.pop_back();//std::remove doesn't remove element but swaps it to the end
                }
                result.swapNeighborhood(i, near, empty, false);
            } else {
                for (index j : near) {
                    if (j >= n) ERROR("Node ", j, " prospective neighbor of ", i, " does not actually exist. Oops.");
                    if(radii[j] > radii[i] || (radii[j] == radii[i] && angles[j] < angles[i]))
                        result.addHalfEdge(i,j);
                }
            }
        }
        threadtimers[id].stop();
    }
    timer.stop();
    INFO("Generating Edges took ", timer.elapsedMilliseconds(), " milliseconds.");
    return result.toGraph(!directSwap, true);
}

Graph HyperbolicGenerator::generate(const vector<double> &angles, const vector<double> &radii, double R, double T) {
    if (T < 0) throw std::runtime_error("Temperature cannot be negative.");
    if (T == 0) return generateCold(angles, radii, R);
    assert(T > 0);

    /**
     * fill Quadtree
     */
    Aux::Timer timer;
    timer.start();
    index n = angles.size();
    assert(radii.size() == n);

    assert(alpha > 0);
    Quadtree<index,false> quad(R, theoreticalSplit, alpha, capacity, balance);

    for (index i = 0; i < n; i++) {
        assert(radii[i] < R);
        quad.addContent(i, angles[i], radii[i]);
    }

    quad.trim();
    timer.stop();
    INFO("Filled Quadtree, took ", timer.elapsedMilliseconds(), " milliseconds.");

    assert(quad.size() == n);

    bool anglesSorted = std::is_sorted(angles.begin(), angles.end());

    //now define lambda
    double beta = 1/T;
    assert(beta == beta);
    auto edgeProb = [beta, R](double distance) -> double {return 1 / (exp(beta*(distance-R)/2)+1);};

    //get Graph
    GraphBuilder result(n, false, false);//no direct swap with probabilistic graphs
    count totalCandidates = 0;
    #pragma omp parallel for reduction(+: totalCandidates)
    for (omp_index i = 0; i < static_cast<omp_index>(n); i++) {
        vector<index> near;
        totalCandidates += quad.getElementsProbabilistically(HyperbolicSpace::polarToCartesian(angles[i], radii[i]), edgeProb, anglesSorted, near);
        for (index j : near) {
            if (j >= n) ERROR("Node ", j, " prospective neighbor of ", i, " does not actually exist. Oops.");
            if (j > i) {
                result.addHalfEdge(i, j);
            }
        }

    }
    DEBUG("Candidates tested: ", totalCandidates);
    return result.toGraph(true, true);

}
} // namespace NetworKit
