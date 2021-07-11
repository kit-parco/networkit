/*
 * SelectiveCommunityDetector.cpp
 *
 *  Created on: 15.05.2013
 *      Author: cls, Yassine Marrakchi
 */

#include <networkit/scd/SelectiveCommunityDetector.hpp>

namespace NetworKit {

SelectiveCommunityDetector::SelectiveCommunityDetector(const Graph &g) : g(&g) {}

std::map<node, std::set<node>> SelectiveCommunityDetector::run(const std::set<node> &seeds) {
    std::map<node, std::set<node>> result;

    for (node s : seeds) {
        result[s] = expandOneCommunity(s);
    }

    return result;
}

std::set<node> SelectiveCommunityDetector::expandOneCommunity(node s) {
    return expandOneCommunity(std::set<node>({s}));
}

} /* namespace NetworKit */
