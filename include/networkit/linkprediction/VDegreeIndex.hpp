// no-networkit-format
/*
 * VDegreeIndex.hpp
 *
 *  Created on: 01.04.2015
 *      Author: Kolja Esders
 */

#ifndef NETWORKIT_LINKPREDICTION_V_DEGREE_INDEX_HPP_
#define NETWORKIT_LINKPREDICTION_V_DEGREE_INDEX_HPP_

#include <networkit/linkprediction/LinkPredictor.hpp>

namespace NetworKit {

/**
 * @ingroup linkprediction
 *
 * Index that simply returns the degree of the second given node.
 */
class VDegreeIndex final : public LinkPredictor {
  /**
   * Returns the degree of the second node provided, namely @a v.
   * @param u First node
   * @param v Second node
   * @return the degree of the second node provided, namely @a v
   */
  double runImpl(node, node v) override {
    return G->degree(v);
  }

public:
  using LinkPredictor::LinkPredictor;

};

} // namespace NetworKit

#endif // NETWORKIT_LINKPREDICTION_V_DEGREE_INDEX_HPP_
