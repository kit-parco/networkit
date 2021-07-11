// no-networkit-format
#include <algorithm>

#include <networkit/community/PartitionIntersection.hpp>

NetworKit::Partition NetworKit::PartitionIntersection::calculate(const Partition &zeta, const NetworKit::Partition &eta) {
    Partition result(std::max(zeta.numberOfElements(), eta.numberOfElements()));
    result.setUpperBound(zeta.upperBound()*eta.upperBound());
    zeta.parallelForEntries([&](node u, index s) {
        if (zeta.contains(u) && eta.contains(u)) {
            result[u] = s * eta.upperBound() + eta[u];
        }
    });
    result.compact();
    return result;
}
