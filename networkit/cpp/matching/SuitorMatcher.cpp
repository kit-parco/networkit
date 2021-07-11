/*
 *  SuitorMatcher.cpp
 *
 *  Created on: 27.08.2019
 *  Authors: Michal Boron     <michal.s.boron@gmail.com>
 *           Eugenio Angriman <angrimae@hu-berlin.de>
 */

#include <algorithm>
#include <atomic>
#include <stdexcept>

#include <networkit/matching/SuitorMatcher.hpp>

namespace NetworKit {

SuitorMatcher::SuitorMatcher(const Graph &G, bool sortSuitor, bool checkSortedEdges)
    : Matcher(G), sortSuitor(sortSuitor) {

    if (G.numberOfSelfLoops() > 0)
        throw std::runtime_error("This algorithm does not graphs with self-loops.");

    if (G.isDirected())
        throw std::runtime_error("This algorithm does not support directed graphs.");

    if (sortSuitor && checkSortedEdges && !hasEdgesSortedByWeight(G))
        throw std::runtime_error("Edges are not sorted by weight");
}

bool SuitorMatcher::hasEdgesSortedByWeight(const Graph &G) {
    std::atomic<bool> isSorted{true};

    G.parallelForNodes([&](node u) {
        if (G.degree(u) < 2 || !isSorted.load(std::memory_order_relaxed))
            return;

        bool uIsSorted;
        if (G.isWeighted())
            uIsSorted = std::is_sorted(
                G.weightNeighborRange(u).begin(), G.weightNeighborRange(u).end(),
                [](const auto &n1, const auto &n2) {
                    return n1.second == n2.second ? n1.first < n2.first : n1.second > n2.second;
                });
        else
            uIsSorted = std::is_sorted(G.neighborRange(u).begin(), G.neighborRange(u).end(),
                                       [](node n1, node n2) { return n1 < n2; });

        if (!uIsSorted)
            isSorted.store(false, std::memory_order_relaxed);
    });

    return isSorted.load(std::memory_order_relaxed);
}

void SuitorMatcher::findSuitor(node current) {
    bool done = false;

    do {

        node partner = suitor[current];
        edgeweight heaviest = ws[current];

        G->forNeighborsOf(current, [&](const node v, const edgeweight weight) {
            if ((weight > heaviest || (weight == heaviest && v < partner))
                && (weight > ws[v] || (weight == ws[v] && current < suitor[v]))) {
                partner = v;
                heaviest = weight;
            }
        });

        done = true;

        if (partner != none
            && (heaviest > ws[partner] || (heaviest == ws[partner] && current < suitor[partner]))) {
            node y = suitor[partner];
            suitor[partner] = current;
            ws[partner] = heaviest;

            // if the current vertex already has a suitor
            if (y != none) {
                current = y;
                done = false;
            }
        }
    } while (!done);
}

void SuitorMatcher::findSortSuitor(node current) {
    bool done = false;

    do {
        node partner = suitor[current];
        edgeweight heaviest = ws[current];

        for (index &iter = neighborIterators[current]; iter < G->degree(current); ++iter) {
            const node v = G->getIthNeighbor(current, iter);
            const edgeweight weight = G->getIthNeighborWeight(current, iter);
            if ((weight > heaviest || (weight == heaviest && v < partner))
                && (weight > ws[v] || (weight == ws[v] && current < suitor[v]))) {
                partner = v;
                heaviest = weight;
                ++iter;
                break;
            }
        }

        done = true;

        if (partner != none
            && (heaviest > ws[partner] || (heaviest == ws[partner] && current < suitor[partner]))) {
            const node y = suitor[partner];
            suitor[partner] = current;
            ws[partner] = heaviest;

            // if the current vertex already has a suitor
            if (y != none) {
                current = y;
                done = false;
            }
        }
    } while (!done);
}

void SuitorMatcher::run() {
    const auto n = G->upperNodeIdBound();
    suitor.assign(n, none);
    ws.assign(n, 0);

    if (sortSuitor)
        neighborIterators.assign(n, 0);

    if (sortSuitor)
        G->forNodes([&](node u) { findSortSuitor(u); });
    else
        G->forNodes([&](node u) { findSuitor(u); });

    G->parallelForNodes([&suitor = suitor, &M = M](node u) {
        if (suitor[u] == none) {
            if (M.isMatched(u))
                M.unmatch(u, M.mate(u));
        } else if (u < suitor[u]) // Ensure we match a pair of nodes only once
            M.match(u, suitor[u]);
    });

    hasRun = true;
}

} // namespace NetworKit
