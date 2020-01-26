/*
 * ParallelPartitionCoarsening.hpp
 *
 *  Created on: 03.07.2014
 *      Author: cls
 */

#ifndef NETWORKIT_COARSENING_PARALLEL_PARTITION_COARSENING_HPP_
#define NETWORKIT_COARSENING_PARALLEL_PARTITION_COARSENING_HPP_

#include <networkit/Globals.hpp>
#include <networkit/coarsening/GraphCoarsening.hpp>
#include <networkit/structures/Partition.hpp>

namespace NetworKit {

/**
 * @ingroup coarsening
 */
class ParallelPartitionCoarsening final : public GraphCoarsening {
public:
    /**
     * All nodes in G have to be assigned a valid subset (!= none) by zeta
     * @param G
     * @param zeta
     * @param parallel
     */
    ParallelPartitionCoarsening(const Graph &G, const Partition &zeta, bool parallel = true);

    void run() override;

private:
    const Partition& zeta;
    bool parallel;
};

} /* namespace NetworKit */

#endif // NETWORKIT_COARSENING_PARALLEL_PARTITION_COARSENING_HPP_
