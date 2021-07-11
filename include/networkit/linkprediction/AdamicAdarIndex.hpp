// no-networkit-format
/*
 * AdamicAdarIndex.hpp
 *
 *  Created on: 25.03.2015
 *      Author: Kolja Esders
 */

#ifndef NETWORKIT_LINKPREDICTION_ADAMIC_ADAR_INDEX_HPP_
#define NETWORKIT_LINKPREDICTION_ADAMIC_ADAR_INDEX_HPP_

#include <networkit/linkprediction/LinkPredictor.hpp>

namespace NetworKit {

/**
 * @ingroup linkprediction
 *
 * Implementation of the Adamic/Adar Index.
 * The index sums up the reciprocals of the logarithm of the degree of all
 * common neighbors of u and v.
 */
class AdamicAdarIndex final : public LinkPredictor {
  /**
   * Returns the Adamic/Adar Index of the given node-pair (@a u, @a v).
   * @param u First node
   * @param v Second node
   * @return the Adamic/Adar Index of the given node-pair (@a u, @a v)
   */
  double runImpl(node u, node v) override;

public:
  using LinkPredictor::LinkPredictor;

};

} // namespace NetworKit

#endif // NETWORKIT_LINKPREDICTION_ADAMIC_ADAR_INDEX_HPP_
