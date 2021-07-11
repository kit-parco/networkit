/*
 * CurveballUniformTradeGenerator.cpp
 *
 *  Created on: 01.06.2019
 *    Author: Manuel Penschuck <networkit@manuel.jetzt>
 */
#include <omp.h>

#include <algorithm>
#include <numeric>
#include <vector>

#include <networkit/Globals.hpp>
#include <networkit/auxiliary/Parallel.hpp>
#include <networkit/auxiliary/Random.hpp>
#include <networkit/auxiliary/SignalHandling.hpp>
#include <networkit/graph/GraphTools.hpp>
#include <networkit/randomization/DegreePreservingShuffle.hpp>

namespace NetworKit {

// The algorithms effectively computes for each degree (in the undirected case it's a scalar, in
// the directed case it's a scalar pair) the set of nodes which have this degree and randomly
// reassigns ids *within* these sets. Since directed and undirected degrees have a different
// signature, the algorithm itself is templated and compares degrees only using the < and the
// == operator. Those implementation details are hidden in the DegreePreservingShuffleDetails
// namespace:
namespace DegreePreservingShuffleDetails {
struct DirectedDegree {
    count in;
    count out;
    bool operator<(const DirectedDegree &o) const {
        return std::tie(in, out) < std::tie(o.in, o.out);
    }
    bool operator==(const DirectedDegree &o) const {
        return std::tie(in, out) == std::tie(o.in, o.out);
    }
};

template <typename DegreeT>
struct NodeDegree {
    node id;
    union {
        DegreeT degree;
        node alternativeId;
    };

    NodeDegree() {}
};

using UndirectedNodeDegree = NodeDegree<count>;
using DirectedNodeDegree = NodeDegree<DirectedDegree>;

template <typename DegreeT>
static std::vector<index> computePermutation(std::vector<NodeDegree<DegreeT>> &nodeDegrees) {
    const auto n = nodeDegrees.size();
    using NodeDegreeType = NodeDegree<DegreeT>;

    // Allocate memory for the permutation
    std::vector<index> permutation(n);
#ifndef NDEBUG
    std::fill_n(permutation.begin(), n, std::numeric_limits<node>::max());
#endif

    // We sort by degree and thereby move nodes of equal degrees into consecutive groups
    Aux::Parallel::sort(
        begin(nodeDegrees), end(nodeDegrees),
        [](const NodeDegreeType a, const NodeDegreeType b) { return a.degree < b.degree; });

    // Now shuffle each group
    count numChanges = 0;
    const auto numThreadsToRequest = std::max(1, std::min<int>(omp_get_max_threads(), n / 50000));
#pragma omp parallel num_threads(numThreadsToRequest) reduction(+ : numChanges)
    {
        const auto tid = omp_get_thread_num();
        const auto numThreads = omp_get_num_threads();

        const auto chunk_size = (n + numThreads - 1) / numThreads;
        auto chunkBegin = nodeDegrees.begin() + std::min<size_t>(chunk_size * tid, n);
        auto chunkEnd = nodeDegrees.begin() + std::min<size_t>(chunk_size * (tid + 1), n);

        using Iterator = decltype(chunkBegin);
        auto firstDifferent = [](Iterator begin, Iterator end, DegreeT degree) {
            return std::find_if_not(begin, end,
                                    [=](NodeDegreeType &x) { return degree == x.degree; });
        };

        if (chunkBegin < chunkEnd) {
            // find real start/end of chunk (i.e., consider groups which extends over a boundary)
            if (chunkBegin != nodeDegrees.begin()) {
                // if a group started in the last chunk continues into this thread's chunk,
                // our first elements belong to the other thread and we have to skip over them.
                chunkBegin = firstDifferent(chunkBegin, chunkEnd, std::prev(chunkBegin)->degree);
            }

            chunkEnd = firstDifferent(chunkEnd, nodeDegrees.end(), std::prev(chunkEnd)->degree);

// we will override the degree values in nodeDegrees. This barrier prevents races between
// those writes and the previous search of chunk limits
#pragma omp barrier

            Aux::SignalHandler handler;

            auto &prng = Aux::Random::getURNG();
            for (auto it = chunkBegin; it != chunkEnd;) {
                handler.assureRunning();

                // find end of our current group
                auto groupEnd = firstDifferent(it + 1, chunkEnd, it->degree);

                // copy id (to make Fischer-Yates work)
                for (auto x = it; x != groupEnd; ++x)
                    x->alternativeId = x->id;

                // go through group and randomly shuffle ids using a Fisher-Yates shuffle
                for (size_t len = std::distance(it, groupEnd) - 1; it != groupEnd; ++it, --len) {
                    const auto rand_idx =
                        len == 0 ? 0 : std::uniform_int_distribution<size_t>(0, len)(prng);

                    permutation[it->id] = it[rand_idx].alternativeId;
                    it[rand_idx].alternativeId = it->alternativeId;

                    numChanges++;
                }
            }
        }
    };

    assert(numChanges == n);
    assert(!std::any_of(permutation.cbegin(), permutation.cend(),
                        [](node i) { return i == std::numeric_limits<node>::max(); }));
    assert(std::accumulate(permutation.cbegin(), permutation.cend(), 0llu)
           == static_cast<size_t>(n - 1) * (n) / 2);

    return permutation;
}
} // namespace DegreePreservingShuffleDetails

DegreePreservingShuffle::DegreePreservingShuffle(const Graph &G) : G(&G) {}
DegreePreservingShuffle::~DegreePreservingShuffle() = default;

std::string DegreePreservingShuffle::toString() const {
    return "DegreePreservingShuffle";
}

void DegreePreservingShuffle::run() {
    const auto n = G->numberOfNodes();

    if (G->isDirected()) {
        // generate sequence of tuple (u, deg(u)) for each u
        std::vector<DegreePreservingShuffleDetails::DirectedNodeDegree> nodeDegrees(n);

        G->parallelForNodes([&](const node u) {
            nodeDegrees[u].id = u;
            nodeDegrees[u].degree = {G->degreeIn(u), G->degreeOut(u)};
        });

        permutation = DegreePreservingShuffleDetails::computePermutation(nodeDegrees);

    } else {
        // generate sequence of tuple (u, deg(u)) for each u
        std::vector<DegreePreservingShuffleDetails::UndirectedNodeDegree> nodeDegrees(n);

        G->parallelForNodes([&](const node u) {
            nodeDegrees[u].id = u;
            nodeDegrees[u].degree = G->degree(u);
        });

        permutation = DegreePreservingShuffleDetails::computePermutation(nodeDegrees);
    }
}

Graph DegreePreservingShuffle::getGraph() const {
    const auto n = G->numberOfNodes();
    assert(permutation.size() == n);
    // localPerm introduced due to a bug in AppleClang 11.03. Direct access to permutation
    // vector creates a seqfault in the compiler.
    const std::vector<node> &localPerm = permutation;

    return GraphTools::getRemappedGraph(*G, n, [&localPerm](node u) { return localPerm[u]; });
}

} // namespace NetworKit
