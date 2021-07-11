// no-networkit-format
#ifndef NETWORKIT_GENERATORS_LFR_GENERATOR_HPP_
#define NETWORKIT_GENERATORS_LFR_GENERATOR_HPP_

#include <networkit/base/Algorithm.hpp>
#include <networkit/generators/StaticGraphGenerator.hpp>
#include <networkit/graph/Graph.hpp>
#include <networkit/structures/Partition.hpp>

#include <tlx/define/deprecated.hpp>

namespace NetworKit {

/**
 * @ingroup generators
 * The LFR clustered graph generator as introduced by Andrea Lancichinetti, Santo Fortunato, and Filippo Radicchi.
 *
 * The community assignment follows the algorithm described in
 * "Benchmark graphs for testing community detection algorithms". The edge generation is however taken from their follow-up publication
 * "Benchmarks for testing community detection algorithms on directed and weighted graphs with overlapping communities". Parts of the
 * implementation follow the choices made in their implementation which is available at https://sites.google.com/site/andrealancichinetti/software
 * but other parts differ, for example some more checks for the realizability of the community and degree size distributions are done
 * instead of heavily modifying the distributions.
 *
 * The edge-switching markov-chain algorithm implementation in NetworKit is used which is different from the implementation in the original LFR benchmark.
 */
class LFRGenerator final : public Algorithm, public StaticGraphGenerator {
public:
    /**
     * Initialize the LFR generator for @a n nodes.
     *
     * @note You need to set a degree sequence, a community size sequence and a mu using the additionally provided set- or generate-methods.
     *
     * @param n The number of nodes.
     */
    LFRGenerator(count n);

    /**
     * Set the given degree sequence.
     *
     * @param degreeSequence The degree sequence that shall be used by the generator
     */
    void setDegreeSequence(std::vector<count> degreeSequence);

    /**
     * Generate and set a power law degree sequence using the given average and maximum degree with the given exponent.
     *
     * @param avgDegree The average degree that shall be reached.
     * @param maxDegree The maximum degree that shall be generated.
     * @param nodeDegreeExp The (negative) exponent of the powerlaw degree sequence.
     */
    void generatePowerlawDegreeSequence(count avgDegree, count maxDegree, double nodeDegreeExp);

    /**
     * Set the given community size sequence.
     *
     * @param communitySizeSequence The community sizes that shall be used.
     */
    void setCommunitySizeSequence(std::vector<count> communitySizeSequence);

    /**
     * Set the partition, this replaces the community size sequence and the random assignment of the nodes to communities.
     *
     * @param zeta The partition to use
     */
    void setPartition(Partition zeta);

    /**
     * Generate a powerlaw community size sequence with the given minimum and maximum size and the given exponent.
     *
     * @param minCommunitySize The minimum community size to generate
     * @param maxCommunitySize The maximum community size to generate
     * @param communitySizeExp The (negative) exponent of the power law community size sequence
     */
    void generatePowerlawCommunitySizeSequence(count minCommunitySize, count maxCommunitySize, double communitySizeExp);

    /**
     * Set the mixing parameter, this is the fraction of neighbors of each node that do not belong to the node's own community.
     *
     * @param mu The mixing parameter that shall be set.
     */
    void setMu(double mu);

    /**
     * Set the mixing parameter separately for each node. This is for each node the fraction of neighbors that do not belong to the node's own community.
     *
     * @param mu The mixing parameter for each node.
     */
    void setMu(const std::vector<double> & mu);

    /**
     * Set the internal degree of each node using a binomial distribution such that the expected mixing parameter is the given @a mu.
     *
     * The mixing parameter is for each node the fraction of neighbors that do not belong to the node's own community.
     *
     * @param mu The expected mu that shall be used.
     */
    void setMuWithBinomialDistribution(double mu);

    /**
     * Generates the graph and the community structure. The algorithm is not parallel.
     */
    void run() override;

    /**
     * Generates and returns the graph.
     *
     * @return The generated graph.
     */
    Graph generate() override;

    /**
     * Returns (a copy of) the generated graph.
     *
     * @return The generated graph.
     */
    Graph getGraph() const;

    /**
     * Returns the generated graph using move semantics.
     *
     * @return The generated graph.
     */
    Graph&& getMoveGraph();

    /**
     * Returns (a copy of) the generated partition
     *
     * @return The generated graph.
     */
    Partition getPartition() const;

    /**
     * Returns the generated partition using move semantics.
     *
     * @return The generated partition.
     */
    Partition&& getMovePartition();

    /**
     * The name and parameters of the generator
     */
    std::string TLX_DEPRECATED(toString() const override);

    /**
     * If the algorithm uses parallelism (no)
     *
     * @return false, only minor parts are parallelized
     */
    bool TLX_DEPRECATED(isParallel() const override);

private:
    /*
     * These methods might be overridden by a sub-class which could use a different model or generator in order to generate the parts of the graph.
     */
    std::vector<std::vector<node>> assignNodesToCommunities();
    Graph generateIntraClusterGraph(std::vector<count> intraDegreeSequence, const std::vector<node> &localToGlobalNode);
    Graph generateInterClusterGraph(const std::vector<count> &externalDegreeSequence);

    count n;
    bool hasDegreeSequence;
    std::vector<count> degreeSequence;
    bool hasCommunitySizeSequence;
    std::vector<count> communitySizeSequence;
    bool hasInternalDegreeSequence;
    std::vector<count> internalDegreeSequence;
    bool hasGraph;
    Graph G;
    bool hasPartition;
    Partition zeta;
};

} // namespace NetworKit

#endif // NETWORKIT_GENERATORS_LFR_GENERATOR_HPP_
