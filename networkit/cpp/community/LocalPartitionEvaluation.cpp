// no-networkit-format
#include <networkit/auxiliary/Log.hpp>
#include <networkit/community/LocalPartitionEvaluation.hpp>

NetworKit::LocalPartitionEvaluation::LocalPartitionEvaluation(const Graph &G, const Partition &P) : G(&G), P(&P) {
    if (P.upperBound() > 2 * G.upperNodeIdBound()) {
        WARN("The upper bound of the partition ", P.upperBound(), " is much higher than the maximum node id: ",
             G.upperNodeIdBound(), ". This might result in high running times and high memory consumption or even crash the whole program.");
    }
}

