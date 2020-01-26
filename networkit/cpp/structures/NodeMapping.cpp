/*
 * NodeMapping.cpp
 *
 * Created: 2019-03-28
 * Author: Armin Wiebigke
 */

#include <networkit/structures/NodeMapping.hpp>

namespace NetworKit {

NodeMapping::NodeMapping(const NetworKit::Graph &G)
        : globalToLocal(G.upperNodeIdBound(), none) {}

NodeMapping::NodeMapping(count globalUpperBound) : globalToLocal(globalUpperBound, none) {}

void NodeMapping::reset() {
    for (node v : localToGlobal) {
        globalToLocal[v] = none;
    }
    localToGlobal.clear();
}

void NodeMapping::reset(index end) {
    for (node v = end; v < localToGlobal.size(); ++v)
        globalToLocal[localToGlobal[v]] = none;
    localToGlobal.resize(end);
}

} /* namespace NetworKit */
