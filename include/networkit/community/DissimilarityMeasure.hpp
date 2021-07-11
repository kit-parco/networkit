// no-networkit-format
/*
 * DissimilarityMeasure.hpp
 *
 *  Created on: 19.01.2013
 *      Author: Christian Staudt
 */

#ifndef NETWORKIT_COMMUNITY_DISSIMILARITY_MEASURE_HPP_
#define NETWORKIT_COMMUNITY_DISSIMILARITY_MEASURE_HPP_

#include <networkit/graph/Graph.hpp>
#include <networkit/structures/Cover.hpp>
#include <networkit/structures/Partition.hpp>

namespace NetworKit {

/**
 * @ingroup community
 * Base class for all clustering dissimilarity measures.
 */
class DissimilarityMeasure {

public:

    virtual double getDissimilarity(const Graph& G, const Partition& first, const Partition& second) = 0;


    virtual double getDissimilarity(const Graph &G, const Cover &first, const Cover &second);
};

} /* namespace NetworKit */
#endif // NETWORKIT_COMMUNITY_DISSIMILARITY_MEASURE_HPP_
