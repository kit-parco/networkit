// no-networkit-format
/*
 * HashingOverlapper.hpp
 *
 *  Created on: 31.01.2013
 *      Author: Christian Staudt
 */

#ifndef NETWORKIT_OVERLAP_HASHING_OVERLAPPER_HPP_
#define NETWORKIT_OVERLAP_HASHING_OVERLAPPER_HPP_

#include <functional>

#include <networkit/overlap/Overlapper.hpp>

namespace NetworKit {

/**
 * @ingroup overlap
 * Determines the overlap of multiple partitions by hashing partition identifiers.
 */
class HashingOverlapper final : public Overlapper {

public:

    Partition run(const Graph& G, const std::vector<Partition>& clusterings) override;

};

} /* namespace NetworKit */
#endif // NETWORKIT_OVERLAP_HASHING_OVERLAPPER_HPP_
