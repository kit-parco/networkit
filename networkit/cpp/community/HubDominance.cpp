// no-networkit-format
/*
 * HubDominance implementation
 *
 * Created: 2014-08-13
 * Author: Michael Hamann
 */

#include <networkit/community/HubDominance.hpp>
#include <networkit/community/PartitionHubDominance.hpp>
#include <networkit/community/CoverHubDominance.hpp>

double NetworKit::HubDominance::getQuality(const NetworKit::Partition &zeta, const NetworKit::Graph &G) {
    PartitionHubDominance phd(G, zeta);
    phd.run();
    return phd.getUnweightedAverage();
}

double NetworKit::HubDominance::getQuality(const NetworKit::Cover &zeta, const NetworKit::Graph &G) {
    CoverHubDominance chd(G, zeta);
    chd.run();
    return chd.getUnweightedAverage();
}
