// no-networkit-format
#ifndef NETWORKIT_COMMUNITY_LOCAL_COMMUNITY_EVALUATION_HPP_
#define NETWORKIT_COMMUNITY_LOCAL_COMMUNITY_EVALUATION_HPP_

#include <vector>

#include <networkit/Globals.hpp>
#include <networkit/base/Algorithm.hpp>

namespace NetworKit {

/**
 * Virtual base class of all evaluation methods for a single clustering which is based on the evaluation of single clusters.
 * This is the base class both for Partitions as well as for Covers.
 */
class LocalCommunityEvaluation : public Algorithm {
public:
    /**
     * Default destructor for the virtual base class
     */
    ~LocalCommunityEvaluation() override = default;

    /**
     * Get the average value weighted by cluster size.
     *
     * @return The weighted average value.
     */
    double getWeightedAverage() const { assureFinished(); return weightedAverage; };

    /**
     * Get the (unweighted) average value of all clusters.
     *
     * @return The unweighted average value.
     */
    double getUnweightedAverage() const { assureFinished(); return unweightedAverage; };

    /**
     * Get the maximum value of all clusters.
     *
     * @return The maximum value.
     */
    double getMaximumValue() const { assureFinished(); return maximumValue; };

    /**
     * Get the minimum value of all clusters.
     *
     * @return The minimum value.
     */
    double getMinimumValue() const { assureFinished(); return minimumValue; };

    /**
     * Get the value of the specified cluster.
     *
     * @param i The cluster to get the value for.
     * @return The value of cluster @a i.
     */
    double getValue(index i) const { assureFinished(); return values[i]; };

    /**
     * Get the values of all clusters.
     *
     * @return The values of all clusters.
     */
    std::vector<double> getValues() const { assureFinished(); return values; };

    /**
     * If small values are better (otherwise large values are better).
     *
     * @return If small values are better.
     */
    virtual bool isSmallBetter() const = 0;
protected:
    double weightedAverage;
    double unweightedAverage;
    double maximumValue;
    double minimumValue;
    std::vector<double> values;
};

}

#endif // NETWORKIT_COMMUNITY_LOCAL_COMMUNITY_EVALUATION_HPP_
