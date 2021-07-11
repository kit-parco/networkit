// no-networkit-format
/*
 * MocnikGenerator.cpp
 *
 * Created on: July 7, 2018
 * Author: Franz-Benjamin Mocnik <mail@mocnik-science.net>
 */

#include <algorithm>
#include <cmath>
#include <map>
#include <utility>
#include <vector>

#include <networkit/auxiliary/Random.hpp>
#include <networkit/generators/MocnikGenerator.hpp>

namespace NetworKit {

MocnikGenerator::MocnikGenerator(count dim, count n, double k, bool weighted): dim(dim), weighted(weighted) {
    ns.push_back(n);
    ks.push_back(k);
}

MocnikGenerator::MocnikGenerator(count dim, std::vector<count> ns, double k, bool weighted)
    : dim(dim), ns(std::move(ns)), weighted(weighted) {
    ks.resize(this->ns.size(), k);
}

MocnikGenerator::MocnikGenerator(count dim, std::vector<count> ns, std::vector<double> ks,
                                 bool weighted)
    : dim(dim), ns(std::move(ns)), ks(std::move(ks)), weighted(weighted) {}

MocnikGenerator::MocnikGenerator(count dim, count n, double k, std::vector<double> weighted)
    : dim(dim), weighted(true), relativeWeights(std::move(weighted)) {
    ns.push_back(n);
    ks.push_back(k);
}

MocnikGenerator::MocnikGenerator(count dim, std::vector<count> ns, double k,
                                 std::vector<double> weighted)
    : dim(dim), ns(std::move(ns)), weighted(true), relativeWeights(std::move(weighted)) {
    ks.resize(this->ns.size(), k);
}

MocnikGenerator::MocnikGenerator(count dim, std::vector<count> ns, std::vector<double> ks,
                                 std::vector<double> weighted)
    : dim(dim), ns(std::move(ns)), ks(std::move(ks)), weighted(true),
      relativeWeights(std::move(weighted)) {}

// TOOLS

/**
 * For a vector of pairs, get the pair with the maximal second component
 */
template<typename KeyType, typename ValueType> std::pair<KeyType,ValueType> getMax(const std::map<KeyType,ValueType> &x) {
    using pairtype=std::pair<KeyType,ValueType>;
    return *std::max_element(x.begin(), x.end(), [](const pairtype &p1, const pairtype &p2) {
        return p1.second < p2.second;
    });
}

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

// LAYER STATE

void MocnikGenerator::initCellArray(MocnikGenerator::LayerState &s, count numberOfCellsPerDimension) {
    s.aMax = numberOfCellsPerDimension;
    for (count j = 0; j < std::pow(s.aMax, dim); j++) {
        NodeCollection tmp;
        s.a.push_back(tmp);
    }
}

MocnikGenerator::NodeCollection MocnikGenerator::getNodes(MocnikGenerator::LayerState &s, int i) {
    return s.a[i];
}

void MocnikGenerator::addNode(MocnikGenerator::LayerState &s, int j) {
    s.a[toIndex(s, nodePositions[j])].push_back(j);
}

int MocnikGenerator::toIndex(MocnikGenerator::LayerState &s, const std::vector<double> &v) {
    std::vector<int> w;
    for (count j = 0; j < v.size(); j++) {
        w.push_back(fmin(floor(v[j] * s.aMax), s.aMax - 1));
    }
    return toIndex(s, w);
}

int MocnikGenerator::toIndex(MocnikGenerator::LayerState &s, const std::vector<int> &v) {
    int x = 0;
    for (count j = v.size() - 1; j < v.size(); j--) {
        x = x * s.aMax + v[j];
    }
    return x;
}

std::vector<int> MocnikGenerator::fromIndex(MocnikGenerator::LayerState &s, int i) {
    std::vector<int> v;
    int i2 = i;
    for (count j = 0; j < dim; j++) {
        int i2New = i2 % s.aMax;
        v.push_back(i2New);
        i2 = (i2 - i2New) / s.aMax;
    }
    return v;
}

std::vector<int> MocnikGenerator::boxSurface(MocnikGenerator::LayerState &s, int i, int r) {
    // test for vanishing r
    if (r == 0) {
        std::vector<int> seResult;
        seResult.push_back(i);
        return seResult;
    }
    // find all boxes
    std::vector<std::vector<int>> se;
    std::vector<int> iV = fromIndex(s, i);
    for (count d = 0; d < dim; d++) {
        std::vector<std::vector<int>> v;
        std::vector<int> tmp;
        v.push_back(tmp);
        for (count j = 0; j < d; j++) {
            std::vector<std::vector<int>> w;
            for (int mu = std::max(iV[j] - r + 1, 0); mu <= std::min(iV[j] + r - 1, s.aMax - 1); mu++) {
                for (std::vector<int> &vElem : v) {
                    std::vector<int> x(vElem);
                    x.push_back(mu);
                    w.push_back(x);
                }
            }
            v = w;
        }
        std::vector<std::vector<int>> w;
        for (std::vector<int> &vElem : v) {
            if (iV[d] - r >= 0) {
                std::vector<int> x(vElem);
                x.push_back(iV[d] - r);
                w.push_back(x);
            }
            if (iV[d] + r < s.aMax) {
                std::vector<int> x(vElem);
                x.push_back(iV[d] + r);
                w.push_back(x);
            }
        }
        v = w;
        for (count j = d + 1; j < dim; j++) {
            std::vector<std::vector<int>> w;
            for (int mu = std::max(iV[j] - r, 0); mu <= std::min(iV[j] + r, s.aMax - 1); mu++) {
                for (std::vector<int> &vElem : v) {
                    std::vector<int> x(vElem);
                    x.push_back(mu);
                    w.push_back(x);
                }
            }
            v = w;
        }
        se.insert(se.end(), v.begin(), v.end());
    }
    // convert the list of grid cells from a multi to a one-dimensional index
    std::vector<int> seResult;
    seResult.reserve(se.size());
    for (std::vector<int> &v : se) {
        seResult.push_back(toIndex(s, v));
    }
    // remove double numbers from the list of grid cells
    sort(seResult.begin(), seResult.end());
    seResult.erase(unique(seResult.begin(), seResult.end()), seResult.end());
    return seResult;
}

std::vector<int> MocnikGenerator::boxVolume(MocnikGenerator::LayerState &s, int i, double r) {
    int r2 = ceil(r * s.aMax);
    std::vector<std::vector<int>> se;
    std::vector<int> tmp;
    se.push_back(tmp);
    std::vector<int> iV = fromIndex(s, i);
    // find all boxes
    for (count d = 0; d < dim; d++) {
        std::vector<std::vector<int>> w;
        for (int mu = std::max(iV[d] - r2, 0); mu <= std::min(iV[d] + r2, s.aMax - 1); mu++) {
            for (std::vector<int> &sElem : se) {
                std::vector<int> x(sElem);
                x.push_back(mu);
                w.push_back(x);
            }
        }
        se = w;
    }
    // convert the list of grid cells from a multi to a one-dimensional index
    std::vector<int> seResult;
    seResult.reserve(se.size());
    for (std::vector<int> &v : se) {
        seResult.push_back(toIndex(s, v));
    }
    return seResult;
}

// EDGE GENERATION

void MocnikGenerator::addEdgesToGraph(Graph &G, count n, double k, double relativeWeight, bool baseLayer) {
    // map vector containing the nodes resp. their positions
    MocnikGenerator::LayerState s;
    initCellArray(s, ceil(std::pow(n / 2, 1./dim) / k));

    // add the nodes to the layer state
    for (count i = 0; i < n; i++) {
        addNode(s, i);
    }

    // create the edges
    count cellMax = std::pow(s.aMax, dim);
    std::vector<std::vector<std::tuple<node, node, double>>> edges(cellMax);
    #pragma omp parallel for
    for (omp_index cell = 0; cell < static_cast<omp_index>(cellMax); cell++) {
        double dmNew, dm;
        edges[cell] = {};
        NodeCollection ns = getNodes(s, cell);
        if (ns.empty()) {
            continue;
        }
        // compute the minimal distance from a node to all other nodes
        std::map<node, double> distMins;
        for (node &i : ns) {
            distMins[i] = -1;
        }
        bool nodesFound = false;
        int r = -1;
        while (!nodesFound || r < getMax(distMins).second * s.aMax) {
            r++;
            for (auto &x : boxSurface(s, cell, r)) {
                NodeCollection ns2 = getNodes(s, x);
                for (node &i : ns) {
                    dm = distMins[i];
                    for (node &j : ns2) {
                        if (i != j) {
                            dmNew = dist(nodePositions[i], nodePositions[j]);
                            if (dmNew < dm || dm == -1) {
                                nodesFound = true;
                                dm = dmNew;
                            }
                        }
                    }
                    distMins[i] = dm;
                }
            }
        }
        // add the edges
        double kdMin, d;
        for (node &i : ns) {
            kdMin = k * distMins[i];
            for (auto &x : boxVolume(s, cell, kdMin)) {
                for (node &j : getNodes(s, x)) {
                    d = dist(nodePositions[i], nodePositions[j]);
                    if (d <= kdMin && i != j) {
                        edges[cell].push_back(std::make_tuple(i, j, d));
                    }
                }
            }
        }
    }

    // add the edges to the graph
    for (count t = 0; t < cellMax; t++) {
        for (auto &e : edges[t]) {
            if (baseLayer || !G.hasEdge(std::get<0>(e), std::get<1>(e))) {
                G.addEdge(std::get<0>(e), std::get<1>(e), std::get<2>(e) * relativeWeight);
            }
        }
    }
}

// GRAPH GENERATION

Graph MocnikGenerator::generate() {
    // create relative weights
    if (relativeWeights.empty()) {
        relativeWeights.resize(ns.size(), 1);
    }

    // assertions
    assert(dim > 0);
    assert(std::all_of(ns.cbegin(), ns.cend(), [] (count n) { return n > 1; }));
    assert(std::all_of(ks.cbegin(), ks.cend(), [] (double k) {return k > 1.0; }));
    assert(!ns.empty());
    assert(ks.size() == ns.size());
    assert(relativeWeights.size() == ns.size());

    // create graph
    Graph G(0, weighted, true);

    // create the nodes
    node curr = 0;
    while (curr < *std::max_element(ns.begin(), ns.end())) {
        std::vector<double> v = {};
        for (count j = 0; j < dim; j++) {
            v.push_back(Aux::Random::real());
        }
        // test whether the new node would be contained in the ball B_{.5}(.5, ..., .5)
        if (norm(v, -.5) <= .5) {
            G.addNode();
            nodePositions.push_back(v);
            curr++;
        }
    }

    // create the edges
    for (count j = 0; j < ns.size(); j++) {
        addEdgesToGraph(G, ns[j], ks[j], relativeWeights[j], j == 0);
    }

    // shrink the graph
    G.shrinkToFit();

    return G;
}

} /* namespace NetworKit */
