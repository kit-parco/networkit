// no-networkit-format
/*
 * PreferentialAttachmentIndex.hpp
 *
 *  Created on: 22.03.2015
 *      Author: Kolja Esders
 */

#ifndef NETWORKIT_LINKPREDICTION_PREFERENTIAL_ATTACHMENT_INDEX_HPP_
#define NETWORKIT_LINKPREDICTION_PREFERENTIAL_ATTACHMENT_INDEX_HPP_

#include <networkit/linkprediction/LinkPredictor.hpp>

namespace NetworKit {

/**
 * @ingroup linkprediction
 *
 * Implementation of the Preferential Attachment Index.
 * The run-method simply calculates the product of the number of nodes in the neighborhoods
 * regarding the given nodes.
 */
class PreferentialAttachmentIndex final : public LinkPredictor {
  /**
   * Returns the product of the cardinalities of the neighborhoods regarding @a u and @a v.
   * @param u First node
   * @param v Second node
   * @return the product of the cardinalities of the neighborhoods regarding @a u and @a v
   */
  double runImpl(node u, node v) override {
    return G->degree(u) * G->degree(v);
  }

public:
  using LinkPredictor::LinkPredictor;

};

} // namespace NetworKit

#endif // NETWORKIT_LINKPREDICTION_PREFERENTIAL_ATTACHMENT_INDEX_HPP_
