// no-networkit-format
/*
 * Volume.cpp
 *
 * Created on: July 7, 2018
 * Author: Franz-Benjamin Mocnik <mail@mocnik-science.net>
 */

#include <networkit/distance/Volume.hpp>
#include <networkit/graph/Graph.hpp>
#include <networkit/graph/GraphTools.hpp>

namespace NetworKit {

std::unordered_map<node, double> Volume::nodesWithinDistance(const Graph &G, double r, node n) {
    std::unordered_map<node, double> ms;
    ms.insert(std::make_pair(n, 0));
    std::vector<node> msToCheck;
    msToCheck.push_back(n);
    double r2;
    while (!msToCheck.empty()) {
        std::vector<node> msToCheckNew;
        for (auto &m : msToCheck) {
            for (node m2 : G.neighborRange(m)) {
                r2 = ms[m] + G.weight(m, m2);
                if (ms.count(m2) == 0) {
                    if (r2 <= r) {
                        ms[m2] = r2;
                        msToCheckNew.push_back(m2);
                    }
                } else {
                    ms[m2] = fmin(ms[m2], r2);
                }
            }
        }
        msToCheck = msToCheckNew;
    }
    return ms;
}

double Volume::volume(const Graph &G, double r, count samples) {
    double x = 0;
    for (count j = 0; j < samples; j++) {
        x += Volume::nodesWithinDistance(G, r, GraphTools::randomNode(G)).size();
    }
    return x / samples;
}

std::vector<double> Volume::volume(const Graph &G, std::vector<double> rs, count samples) {
    std::vector<double> xs(rs.size(), 0);
    double rmax = *std::max_element(std::begin(rs), std::end(rs));
    for (count j = 0; j < samples; j++) {
        std::unordered_map<node, double> ms = Volume::nodesWithinDistance(G, rmax, GraphTools::randomNode(G));
        count i = 0;
        for (const auto &r : rs) {
            for (const auto &it : ms) {
                if (it.second <= r) {
                    xs[i] += 1;
                }
            }
            i++;
        }
    }
    std::vector<double> ys;
    ys.reserve(xs.size());
    for (const auto &x : xs) {
        ys.push_back(x / samples);
    }
    return ys;
}

} /* namespace NetworKit */
