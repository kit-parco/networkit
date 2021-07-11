// no-networkit-format
/*
 * SampledGraphStructuralRandMeasure.hpp
 *
 *  Created on: 01.07.2013
 *      Author: cls
 */

#ifndef NETWORKIT_COMMUNITY_SAMPLED_GRAPH_STRUCTURAL_RAND_MEASURE_HPP_
#define NETWORKIT_COMMUNITY_SAMPLED_GRAPH_STRUCTURAL_RAND_MEASURE_HPP_

#include <networkit/community/DissimilarityMeasure.hpp>

namespace NetworKit {

/**
 * @ingroup community
 * The graph-structural Rand measure assigns a similarity value in [0,1]
 * to two partitions of a graph, by considering connected pairs of nodes.
 * This implementation approximates the index by sampling.
 */
class SampledGraphStructuralRandMeasure final : public DissimilarityMeasure {

public:

    /**
     * Constructs the SampledGraphStructuralRandMeasure. A maximum of @a maxSamples samples are drawn.
     *
     * @param maxSamples The amount of samples to draw.
     */
    SampledGraphStructuralRandMeasure(count maxSamples);

    double getDissimilarity(const Graph& G, const Partition& first, const Partition& second)override;

private:

    count maxSamples;

};

} /* namespace NetworKit */
#endif // NETWORKIT_COMMUNITY_SAMPLED_GRAPH_STRUCTURAL_RAND_MEASURE_HPP_
