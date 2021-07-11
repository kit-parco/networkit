// no-networkit-format
/*
 * ChungLu.hpp
 *
 *  Created on: Dec 23, 2013
 *      Author: Henning
 */

#ifndef NETWORKIT_GENERATORS_CHUNG_LU_GENERATOR_HPP_
#define NETWORKIT_GENERATORS_CHUNG_LU_GENERATOR_HPP_

#include <networkit/generators/StaticDegreeSequenceGenerator.hpp>

namespace NetworKit {

/**
 * @ingroup generators
 * Given an arbitrary degree sequence, the Chung-Lu generative model
 * will produce a random graph with the same expected degree sequence.
 *
 * see Chung, Lu: The average distances in random graphs with given expected degrees
 * and Chung, Lu: Connected Components in Random Graphs with Given Expected Degree Sequences.
 * Aiello, Chung, Lu: A Random Graph Model for Massive Graphs describes a different generative model
 * which is basically asymptotically equivalent but produces multi-graphs.
 *
 * This follows the implementation of Joel Miller and Aric Hagberg's
 * "Efficient Generation of Networks with Given Expected Degrees" (2011)
 * http://aric.hagberg.org/papers/miller-2011-efficient.pdf .
 * It gives a complexity of O(n+m) as opposed to quadratic.
 *
 */

class ChungLuGenerator final: public StaticDegreeSequenceGenerator {
    count sum_deg;
    count n;

public:

    ChungLuGenerator(const std::vector<count>& degreeSequence);

    /**
     * Generates graph with expected degree sequence seq.
     */
    Graph generate() override;
};

} /* namespace NetworKit */
#endif // NETWORKIT_GENERATORS_CHUNG_LU_GENERATOR_HPP_
