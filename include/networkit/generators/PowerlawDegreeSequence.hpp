// no-networkit-format
#ifndef NETWORKIT_GENERATORS_POWERLAW_DEGREE_SEQUENCE_HPP_
#define NETWORKIT_GENERATORS_POWERLAW_DEGREE_SEQUENCE_HPP_

#include <vector>

#include <networkit/Globals.hpp>
#include <networkit/base/Algorithm.hpp>
#include <networkit/graph/Graph.hpp>

namespace NetworKit {

class PowerlawDegreeSequence final : public Algorithm {
public:
    /**
     * Generates a powerlaw degree sequence with the given minimum and maximum degree, the powerlaw exponent gamma.
     *
     * @param minDeg   The minium degree
     * @param maxDeg   The maximum degree
     * @param gamma    The powerlaw exponent
     */
    PowerlawDegreeSequence(count minDeg, count maxDeg, double gamma);

    /**
     * Generates a powerlaw degree sequence that fits the given degree sequence in terms of minimum, maximum and average degree.
     *
     * @param degreeSequence The degree sequence to fit.
     */
    PowerlawDegreeSequence(const std::vector<double>& degreeSequence);

    /**
     * Generates a powerlaw degree sequence that fits the degree sequence of the given graph in terms of minimum, maximum and average degree.
     *
     * @param g The input graph to fit
     */
    PowerlawDegreeSequence(const Graph& g);

    /**
     * Tries to set the minimum degree such that the specified average degree is expected.
     *
     * @throws std::runtime_error If it is not possible to find a minimum degree such that the expected average is @a avgDeg.
     * @param avgDeg   The average degree
     */
    void setMinimumFromAverageDegree(double avgDeg);

    /**
     * Tries to set the powerlaw exponent gamma such that the specified average degree is expected.
     *
     * @param avgDeg The average degree
     * @param minGamma The minimum gamma to use, default: -1
     * @param maxGamma The maximum gamma to use, default: -6
     */
    void setGammaFromAverageDegree(double avgDeg, double minGamma = -1, double maxGamma = -6);

    /**
     * Sets the minimum degree.
     *
     * @param minDeg The degree that shall be set as minimum degree
     */
    void setMinimumDegree(count minDeg);

    /**
     * Gets the minimum degree.
     *
     * @return The minimum degree.
     */
    count getMinimumDegree() const;

    /**
     * Gets the maximum degree.
     *
     * @return The maximum degree.
     */
     count getMaximumDegree() const { return maxDeg; };

     /**
      * Sets the exponent gamma.
      *
      * @param gamma The exponent, must be negative.
      */
    void setGamma(double gamma);

     /**
      * Gets the exponent gamma.
      *
      * @return gamma
      */
    double getGamma() const { return gamma; };

    /**
     * Execute the generation process
     */
    void run() override;

    /**
     * Returns the expected average degree.
     *
     * @return The expected average degree.
     */
    double getExpectedAverageDegree() const;

    /**
     * Returns a degree sequence of @a numNodes degrees with even degree sum.
     *
     * @param numNodes The number of nodes
     * @return The generated degree sequence.
     */
    std::vector<count> getDegreeSequence(count numNodes) const;

    /**
     * Returns a degree drawn at random with a power law distribution
     *
     * @return A degree that follows the generated distribution.
     */
    count getDegree() const;

private:
    count minDeg, maxDeg;
    double gamma;
    std::vector<double> cumulativeProbability;
};

} // namespace NetworKit

#endif // NETWORKIT_GENERATORS_POWERLAW_DEGREE_SEQUENCE_HPP_
