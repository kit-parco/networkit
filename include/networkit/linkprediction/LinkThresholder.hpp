/*
 * LinkThresholder.hpp
 *
 *  Created on: 01.04.2015
 *      Author: Kolja Esders
 */

#ifndef NETWORKIT_LINKPREDICTION_LINK_THRESHOLDER_HPP_
#define NETWORKIT_LINKPREDICTION_LINK_THRESHOLDER_HPP_

#include <networkit/linkprediction/LinkPredictor.hpp>

namespace NetworKit {

/**
 * @ingroup linkprediction
 *
 * Filters given predictions based on some criterion and returns a vector of
 * node-pairs that fulfill the given criterion.
 * This can be used to determine which node-pairs should actually be interpreted
 * as future links and which shouldn't.
 */
namespace LinkThresholder {

/**
 * Returns the node-pairs whose scores are at least equal to the given @a minScore.
 * @param predictions Predictions to filter
 * @param minScore Minimal score that the returned node-pairs should have
 * @return a vector of node-pairs whose scores are at least equal to the given @a minScore
 */
std::vector<std::pair<node, node>>
byScore(const std::vector<LinkPredictor::prediction> &predictions, double minScore);

/**
 * Returns the first @a numLinks highest scored node-pairs.
 * @param predictions Predictions to filter
 * @param numLinks Number of top-scored node-pairs to return
 * @return the first @a numLinks highest scored node-pairs
 */
std::vector<std::pair<node, node>>
byCount(const std::vector<LinkPredictor::prediction> &predictions, count numLinks);

/**
 * Returns the first @a percentageLinks percent of the highest scores node-pairs.
 * @param predictions Predictions to filter
 * @param percentageLinks Percentage of highest scored node-pairs to return
 * @return the first @a percentageLinks percent of the highest scores node-pairs
 */
std::vector<std::pair<node, node>>
byPercentage(const std::vector<LinkPredictor::prediction> &predictions, double percentageLinks);

} // namespace LinkThresholder

} // namespace NetworKit

#endif // NETWORKIT_LINKPREDICTION_LINK_THRESHOLDER_HPP_
