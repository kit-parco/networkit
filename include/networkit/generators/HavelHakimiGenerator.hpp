// no-networkit-format
/*
 * HavelHakimiGenerator.hpp
 *
 *  Created on: Dec 10, 2013
 *      Author: Henning
 */

#ifndef NETWORKIT_GENERATORS_HAVEL_HAKIMI_GENERATOR_HPP_
#define NETWORKIT_GENERATORS_HAVEL_HAKIMI_GENERATOR_HPP_

#include <vector>

#include <networkit/generators/StaticDegreeSequenceGenerator.hpp>
#include <networkit/graph/Graph.hpp>

namespace NetworKit {

/**
 * @ingroup generators
 * Havel-Hakimi algorithm for generating a graph according to a given degree sequence.
 * The sequence, if it is realizable, is reconstructed exactly.
 * The resulting graph usually has a high clustering coefficient.
 * Construction runs in linear time O(m). However, the test if a sequence is realizable
 * is quadratic in the sequence length.
 */
class HavelHakimiGenerator final : public StaticDegreeSequenceGenerator  {
public:
    /**
     * @param[in] sequence Degree sequence to realize. Must be non-increasing.
     * @param[in] ignoreIfRealizable If true, generate the graph even if the degree sequence is not realizable. Some nodes may get lower degrees than requested in the sequence.
     */
    HavelHakimiGenerator(const std::vector<count>& sequence, bool ignoreIfRealizable = false);

    /**
     * Generates degree sequence seq (if it is realizable).
     * @throws std::runtime_error If the sequence is not realizable and ignoreIfRealizable is false.
     * @return Graph with degree sequence seq or modified sequence if ignoreIfRealizable is true and the sequence is not realizable.
     */
    Graph generate() override;
private:
    bool ignoreIfRealizable;
};


} /* namespace NetworKit */
#endif // NETWORKIT_GENERATORS_HAVEL_HAKIMI_GENERATOR_HPP_
