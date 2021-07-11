// no-networkit-format
/*
 * HarmonicCloseness.h
 *
 * Created on: 24.02.2018
 * 		 Author: Eugenio Angriman
 */

#ifndef NETWORKIT_CENTRALITY_HARMONIC_CLOSENESS_HPP_
#define NETWORKIT_CENTRALITY_HARMONIC_CLOSENESS_HPP_

#include <networkit/centrality/Centrality.hpp>

namespace NetworKit {

/**
 * @ingroup centrality
 */
class HarmonicCloseness : public Centrality {
public:
  /**
   * Constructs the HarmonicCloseness class for the given Graph @a G. If
   * the closeness scores should be normalized, then set @a normalized to
   * <code>true</code>. The run() method takes O(nm) time, where n is the number
   * of nodes and m is the number of edges of the graph.
   *
   * @param G The graph.
   * @param normalized Set this parameter to <code>false</code> if scores should
   * not be normalized into an interval of [0, 1]. Normalization only for
   * unweighted graphs.
   *
   */
  HarmonicCloseness(const Graph &G, bool normalized = true);

  /**
   * Computes the harmonic closeness centrality on the graph passed in
   * constructor.
   */
  void run() override;

  /*
   * Returns the maximum possible harmonic closeness centrality that a node can
   * have in a star graph with the same amount of nodes.
   */
  double maximum() override;
};
} // namespace NetworKit

#endif // NETWORKIT_CENTRALITY_HARMONIC_CLOSENESS_HPP_
