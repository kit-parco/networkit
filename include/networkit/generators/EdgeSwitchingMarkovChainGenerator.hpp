#ifndef NETWORKIT_GENERATORS_EDGE_SWITCHING_MARKOV_CHAIN_GENERATOR_HPP_
#define NETWORKIT_GENERATORS_EDGE_SWITCHING_MARKOV_CHAIN_GENERATOR_HPP_

#include <networkit/generators/StaticDegreeSequenceGenerator.hpp>

namespace NetworKit {

/**
 * @ingroup generators
 * Graph generator for generating a random simple graph with exactly the given degree sequence based
 * on the Edge-Switching Markov-Chain method.
 *
 * This implementation is based on the paper
 * "Random generation of large connected simple graphs with prescribed degree distribution" by
 * Fabien Viger and Matthieu Latapy, available at http://www-rp.lip6.fr/~latapy/FV/generation.html,
 * however without preserving connectivity (this could later be added as optional feature).
 *
 * The Havel-Hakami generator is used for the initial graph generation, then the Markov-Chain
 * Monte-Carlo algorithm as described and implemented by Fabien Viger and Matthieu Latapy but
 * without the steps for ensuring connectivity is executed. This should lead to a graph that is
 * drawn uniformly at random from all graphs with the given degree sequence.
 *
 * Note that at most 10 times the number of edges edge swaps are performed (same number as in the
 * aforementioned implementation) and in order to limit the running time, at most 200 times as many
 * attempts to perform an edge swap are made (as certain degree distributions do not allow edge
 * swaps at all).
 */
class EdgeSwitchingMarkovChainGenerator final : public StaticDegreeSequenceGenerator {
public:
    /**
     * Initializes the configuration model generator with the given degree sequence @a sequence.
     *
     * If the sequence cannot be realized, optionally if @a ignoreIfRealizable is true, a graph with
     * a different degree sequence is generated where certain nodes might not have their full
     * degree.
     *
     * @param sequence The wanted degree sequence.
     * @param ignoreIfNotRealizable Ignore if the sequence is not realizable and generate a graph
     *        anyway.
     * @param numSwitchesPerEdge Average number of edge switches (Markov chain steps) per edge
     */
    EdgeSwitchingMarkovChainGenerator(const std::vector<count> &sequence,
                                      bool ignoreIfNotRealizable = false,
                                      count numSwitchesPerEdge = 10);

    /**
     * Generate a graph according to the configuration model.
     *
     * Issues a INFO log message if the wanted number of edge swaps cannot be performed because of
     * the limit of attempts (see in the description of the class for details).
     *
     * @throws std::runtime_error If the sequence cannot be generated and ignoreIfRealizable is
     * false.
     * @return The generated graph
     */
    Graph generate() override;

private:
    bool ignoreIfNotRealizable;
    count numSwitchesPerEdge;
};

} // namespace NetworKit

#endif // NETWORKIT_GENERATORS_EDGE_SWITCHING_MARKOV_CHAIN_GENERATOR_HPP_
