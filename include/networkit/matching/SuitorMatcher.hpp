/*
 *  SuitorMatcher.hpp
 *
 *  Created on: 27.08.2019
 *  Authors: Michal Boron     <michal.s.boron@gmail.com>
 *           Eugenio Angriman <angrimae@hu-berlin.de>
 */

#ifndef NETWORKIT_MATCHING_SUITOR_MATCHER_HPP_
#define NETWORKIT_MATCHING_SUITOR_MATCHER_HPP_

#include <vector>

#include <networkit/graph/Graph.hpp>
#include <networkit/matching/Matcher.hpp>

namespace NetworKit {

/**
 * @ingroup matching
 * Suitor matching finding algorithm.
 */
class SuitorMatcher final : public Matcher {

    static bool hasEdgesSortedByWeight(const Graph &G);

    void findSuitor(node current);
    void findSortSuitor(node current);

public:
    /**
     * Computes a 1/2-approximation of the maximum (weighted) matching of an undirected graph using
     * the Suitor algorithm from Manne and Halappanavar presented in "New Effective Multithreaded
     * Matching Algorithms", IPDPS 2014. The algorithm has two versions: SortSuitor (faster, but
     * only works on graphs with adjacency lists sorted by non-increasing edge weight) and Suitor
     * (works on generic graphs). If using SortSuitor, use GrapTools::sortEdgesByWeight(G, true) to
     * sort the adjacency lists by non-increasing edge weight.
     *
     * @param G An undirected graph.
     * @param sortSuitor If true uses the SortSuitor version, otherwise it uses Suitor.
     * @param checkSortedEdges If true and sortSuitor is true it checks whether the adjacency lists
     * of the input graph are sorted by non-increasing edge weight. If they are not, it throws a
     * std::runtime_error.
     */
    SuitorMatcher(const Graph &G, bool sortSuitor = true, bool checkSortedEdges = false);

    ~SuitorMatcher() override = default;

    /**
     * Runs the algorithm.
     */
    void run() override;

private:
    bool sortSuitor;
    std::vector<node> suitor;
    std::vector<edgeweight> ws;
    std::vector<index> neighborIterators;
};
} // namespace NetworKit

#endif // NETWORKIT_MATCHING_SUITOR_MATCHER_HPP_
