/*
 *  GroupClosenessGrowShrink.cpp
 *
 *  Created on: 19.12.2019
 *      Author: Eugenio Angriman <angrimae@hu-berlin.de>
 */

#include <networkit/centrality/GroupClosenessGrowShrink.hpp>

#include "GroupClosenessGrowShrinkImpl.hpp"

namespace NetworKit {

GroupClosenessGrowShrink::GroupClosenessGrowShrink(const Graph &G, const std::vector<node> &group,
                                                   bool extended, count insertions,
                                                   count maxIterations)
    : G(&G),
      weightedImpl{G.isWeighted() ? std::make_unique<
                       GroupClosenessGrowShrinkDetails::GroupClosenessGrowShrinkImpl<edgeweight>>(
                       G, group, extended, insertions, maxIterations)
                                  : nullptr},
      unweightedImpl{G.isWeighted()
                         ? nullptr
                         : std::make_unique<
                             GroupClosenessGrowShrinkDetails::GroupClosenessGrowShrinkImpl<count>>(
                             G, group, extended, insertions, maxIterations)} {}

GroupClosenessGrowShrink::~GroupClosenessGrowShrink() = default;

void GroupClosenessGrowShrink::run() {
    if (G->isWeighted())
        weightedImpl->run();
    else
        unweightedImpl->run();
}

std::vector<node> GroupClosenessGrowShrink::groupMaxCloseness() const {
    if (G->isWeighted())
        return weightedImpl->groupMaxCloseness();
    else
        return unweightedImpl->groupMaxCloseness();
}

count GroupClosenessGrowShrink::numberOfIterations() const {
    if (G->isWeighted())
        return weightedImpl->numberOfIterations();
    else
        return unweightedImpl->numberOfIterations();
}

} // namespace NetworKit
