/*
 * LouvainMapEquation.hpp
 *
 * Created on: 2019-01-28
 * Author: Armin Wiebigke
 *         Michael Hamann
 *         Lars Gottesbüren
 */

#ifndef NETWORKIT_COMMUNITY_LOUVAIN_MAP_EQUATION_HPP_
#define NETWORKIT_COMMUNITY_LOUVAIN_MAP_EQUATION_HPP_

#include <networkit/auxiliary/SparseVector.hpp>
#include <networkit/auxiliary/SpinLock.hpp>
#include <networkit/community/CommunityDetectionAlgorithm.hpp>

namespace NetworKit {

class LouvainMapEquation : public CommunityDetectionAlgorithm {
private:
    enum class ParallelizationType : uint8_t { None, RelaxMap, Synchronous };

    /**
     * @param[in] G input graph
     * @param[in] hierarchical use recursive coarsening
     * @param[in] maxIterations maximum number of iterations for move phase
     * @param[in] parallelization strategy (default synchronous):
     *            none
     *            relaxmap    -> one lock per community to update cuts
     *            synchronous -> work on stale cuts and volumes, update in second step
     *
     */
    explicit LouvainMapEquation(const Graph &graph, bool hierarchical, count maxIterations,
                                ParallelizationType parallelizationType);

public:
    /**
     * @param[in] G input graph
     * @param[in] hierarchical use recursive coarsening
     * @param[in] maxIterations maximum number of iterations for move phase
     * @param[in] parallelization strategy (default synchronous):
     *            none
     *            relaxmap    -> one lock per community to update cuts
     *            synchronous -> work on stale cuts and volumes, update in second step
     *
     */
    explicit LouvainMapEquation(const Graph &graph, bool hierarchical = false,
                                count maxIterations = 32,
                                const std::string &parallelizationStrategy = "relaxmap")
        : LouvainMapEquation(graph, hierarchical, maxIterations,
                             convertStringToParallelizationType(parallelizationStrategy)) {}

    void run() override;

    std::string toString() const override;

private:
    struct Move {
        node movedNode = none;
        double volume = 0.0;
        index originCluster = none, targetCluster = none;
        double cutUpdateToOriginCluster = 0.0, cutUpdateToTargetCluster = 0.0;

        Move(node n, double vol, index cc, index tc, double cuptoc, double cupttc)
            : movedNode(n), volume(vol), originCluster(cc), targetCluster(tc),
              cutUpdateToOriginCluster(cuptoc), cutUpdateToTargetCluster(cupttc) {}
    };

    static_assert(std::is_trivially_destructible<Move>::value,
                  "LouvainMapEquation::Move struct is not trivially destructible");

    const bool parallel;
    ParallelizationType parallelizationType;

    bool hierarchical;
    count maxIterations;

    std::vector<double> clusterCut, clusterVolume;
    double totalCut, totalVolume;

    // for RelaxMap
    std::vector<Aux::Spinlock> locks;

    // for SLM
    Partition nextPartition;
    std::vector<SparseVector<double>> ets_neighborClusterWeights;

    count localMoving(std::vector<node> &nodes, count iteration);

    count synchronousLocalMoving(std::vector<node> &nodes, count iteration);

    bool tryLocalMove(node u, SparseVector<double> &neighborClusterWeights,
                      std::vector<Move> &moves, bool synchronous);

    bool performMove(node u, double degree, double loopWeight, node currentCluster,
                     node targetCluster, double weightToTarget, double weightToCurrent);

    void aggregateAndApplyCutAndVolumeUpdates(std::vector<Move> &moves);

    void calculateInitialClusterCutAndVolume();

    void runHierarchical();

    /**
     * Calculate the change in the map equation if the node is moved from its current cluster to the
     * target cluster. To simplify the calculation, we remove terms that are constant for all target
     * clusters. As a result, "moving" the node to its current cluster gives a value != 0, although
     * the complete map equation would not change.
     */
    double fitnessChange(node, double degree, double loopWeight, node currentCluster,
                         node targetCluster, double weightToTarget, double weightToCurrent,
                         double totalCutCurrently);

    count idivCeil(count a, count b) const { return (a + b - 1) / b; }

    std::pair<count, count> chunkBounds(count i, count n, count chunkSize) const {
        return std::make_pair(i * chunkSize, std::min(n, (i + 1) * chunkSize));
    }

    ParallelizationType
    convertStringToParallelizationType(const std::string &parallelizationStrategy) const {
        if (parallelizationStrategy == "none")
            return ParallelizationType::None;
        else if (parallelizationStrategy == "relaxmap")
            return ParallelizationType::RelaxMap;
        else if (parallelizationStrategy == "synchronous")
            return ParallelizationType::Synchronous;
        else
            throw std::runtime_error("Invalid parallelization type for map equation Louvain: "
                                     + parallelizationStrategy);
    }

#ifndef NDEBUG
    long double sumPLogPwAlpha = 0;
    long double sumPLogPClusterCut = 0;
    long double sumPLogPClusterCutPlusVol = 0;
    double plogpRel(double w);
    void updatePLogPSums();
    long double mapEquation();
    void checkUpdatedCutsAndVolumesAgainstRecomputation();
#endif // NDEBUG
};

} // namespace NetworKit

#endif // NETWORKIT_COMMUNITY_LOUVAIN_MAP_EQUATION_HPP_
