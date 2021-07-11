// no-networkit-format
/*
 * SampledNodeStructuralRandMeasure.hpp
 *
 *  Created on: 01.07.2013
 *      Author: cls
 */

#ifndef NETWORKIT_COMMUNITY_SAMPLED_NODE_STRUCTURAL_RAND_MEASURE_HPP_
#define NETWORKIT_COMMUNITY_SAMPLED_NODE_STRUCTURAL_RAND_MEASURE_HPP_

#include <networkit/community/DissimilarityMeasure.hpp>

namespace NetworKit {

/**
 * @ingroup community
 * The node-structural Rand measure assigns a similarity value in [0,1]
 * to two partitions of a graph, by considering pairs of nodes.
 * This implementation approximates the index by sampling.
 */
class SampledNodeStructuralRandMeasure final : public DissimilarityMeasure {

public:

    /**
     * Constructs the SampledNodeStructuralRandMeasure. A maximum of @a maxSamples samples are drawn.
     *
     * @param maxSamples The amount of samples to draw.
     */
    SampledNodeStructuralRandMeasure(count maxSamples);

    double getDissimilarity(const Graph& G, const Partition& first, const Partition& second) override;

private:

    count maxSamples;

};

} /* namespace NetworKit */
#endif // NETWORKIT_COMMUNITY_SAMPLED_NODE_STRUCTURAL_RAND_MEASURE_HPP_
