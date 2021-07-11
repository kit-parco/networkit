// no-networkit-format
/*
 * Closeness.hpp
 *
 *  Created on: 03.10.2014
 *     Authors: nemes,
 *              Eugenio Angriman <angrimae@hu-berlin.de>
 */

#ifndef NETWORKIT_CENTRALITY_CLOSENESS_HPP_
#define NETWORKIT_CENTRALITY_CLOSENESS_HPP_

#include <tlx/container/d_ary_addressable_int_heap.hpp>

#include <networkit/centrality/Centrality.hpp>

namespace NetworKit {

enum ClosenessVariant { standard = 0, generalized = 1 };

/**
 * @ingroup centrality
 */
class Closeness : public Centrality {

  public:
    /**
     * Constructs the Closeness class for the given Graph @a G. If the closeness
     * scores should be normalized, then set @a normalized to <code>true</code>.
     * The run() method takes O(nm) time, where n is the number of nodes and m
     * is the number of edges of the graph. NOTICE: the graph has to be
     * connected.
     *
     * @param G The graph.
     * @param normalized Set this parameter to <code>false</code> if scores
     * should not be normalized into an interval of [0, 1]. Normalization only
     * for unweighted graphs.
     *
     */
      Closeness(const Graph &G, bool normalized,
                ClosenessVariant variant = ClosenessVariant::standard);

      /**
       * Old constructor, we keep it for backward compatibility. It computes the
       * standard variant of the closenes.
       *
       * @param G The graph.
       * @param normalized Set this parameter to <code>false</code> if scores
       * should not be normalized into an interval of [0, 1]. Normalization only
       * for unweighted graphs.
       * @param checkConnectedness turn this off if you know the graph is
       * connected.
       *
       */
      Closeness(const Graph &G, bool normalized = true, bool checkConnectedness = true);

      /**
       * Computes closeness cetrality on the graph passed in constructor.
       */
      void run() override;

      /*
       * Returns the maximum possible Closeness a node can have in a graph with
       * the same amount of nodes (=a star)
       */
      double maximum() override { return normalized ? 1. : (1. / (static_cast<double>(G.upperNodeIdBound() - 1))); }

  private:
    ClosenessVariant variant;
    std::vector<std::vector<count>> uDist;
    std::vector<std::vector<double>> dDist;
    std::vector<std::vector<uint8_t>> visited;
    std::vector<uint8_t> ts;

    void checkConnectedComponents() const;
    void bfs();
    void dijkstra();
    void incTS();
    void updateScoreData(node u, count reached, double sum) {
        if (sum > 0) {
            if (variant == ClosenessVariant::standard) {
                scoreData[u] = 1.0/ sum;
            } else {
                scoreData[u] = static_cast<double>(reached - 1) / sum / static_cast<double>(G.numberOfNodes() - 1);
            }
        } else {
            scoreData[u] = 0.;
        }
        if (normalized)
            scoreData[u] *=
                (variant == ClosenessVariant::standard ? G.numberOfNodes()
                                                       : reached) -
                1.;
    }

    struct Compare {
      public:
        Compare(const std::vector<double> &dist_) : dist(dist_) {}

        bool operator()(node x, node y) const { return dist[x] < dist[y]; }

      private:
        const std::vector<double> &dist;
    };

    std::vector<tlx::d_ary_addressable_int_heap<node, 2, Compare>> heaps;
};

} /* namespace NetworKit */

#endif // NETWORKIT_CENTRALITY_CLOSENESS_HPP_
