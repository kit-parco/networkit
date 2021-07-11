// no-networkit-format
/*
 * MocnikGeneratorBasic.cpp
 *
 * Created on: July 7, 2018
 * Author: Franz-Benjamin Mocnik <mail@mocnik-science.net>
 */

#include <networkit/auxiliary/Random.hpp>
#include <networkit/generators/MocnikGeneratorBasic.hpp>

namespace NetworKit {

MocnikGeneratorBasic::MocnikGeneratorBasic(count dim, count n, double k): dim(dim), n(n), k(k) {}

// GEOMETRY

/**
 * Norm of a vector.  The shift applies to every coordinate.
 */
static inline double norm(std::vector<double> &v, double shift) {
    double x = 0;
    for (count j = 0; j < v.size(); j++) {
        x += (v[j] + shift) * (v[j] + shift);
    }
    return sqrt(x);
}

/**
 * Euclidean distance between two vectors
 */
static inline double dist(std::vector<double> &v, std::vector<double> &w) {
    double x = 0;
    for (count j = 0; j < v.size(); j++) {
        x += std::pow(v[j] - w[j], 2);
    }
    x = sqrt(x);
    return x;
}

// GRAPH GENERATION

Graph MocnikGeneratorBasic::generate() {
    // assertions
    assert (dim);
    assert (n);
    assert (k > 1);

    // create graph
    Graph G(0, false, true);

    // create the nodes
    nodePositions.resize(n);
    node curr = 0;
    while (curr < n) {
        std::vector<double> v(dim);
        for (count j = 0; j < dim; j++) {
            v[j] = Aux::Random::real();
        }
        // test whether the new node would be contained in the ball B_{.5}(.5, ..., .5)
        if (norm(v, -.5) < .5) {
            G.addNode();
            nodePositions[curr] = v;
            curr++;
        }
    }

    // create the edges
    double x;
    for (count i = 0; i < n; i++) {
        // compute the minimal distance from a node to all other nodes
        double distMin = -1;
        for (count j = 0; j < n; j++) {
            if (i != j) {
                x = dist(nodePositions[i], nodePositions[j]);
                if (x < distMin || distMin <= -.5) {
                    distMin = x;
                }
            }
        }
        // add the edges
        for (count j = 0; j < n; j++) {
            if (i != j && dist(nodePositions[i], nodePositions[j]) <= k * distMin) {
                G.addEdge(i, j);
            }
        }
    }

    // shrink the graph
    G.shrinkToFit();

    return G;
}

} /* namespace NetworKit */
