// no-networkit-format
/*
 * Overlapper.hpp
 *
 *  Created on: 30.10.2012
 *      Author: Christian Staudt
 */

#ifndef NETWORKIT_OVERLAP_OVERLAPPER_HPP_
#define NETWORKIT_OVERLAP_OVERLAPPER_HPP_

#include <networkit/graph/Graph.hpp>
#include <networkit/structures/Partition.hpp>

namespace NetworKit {

/**
 * @ingroup overlap
 * Abstract base class for algorithms which determine the overlap of multiple partitions.
 */
class Overlapper {

public:

    virtual Partition run(const Graph& G, const  std::vector<Partition>& clusterings) = 0;

};

} /* namespace NetworKit */
#endif // NETWORKIT_OVERLAP_OVERLAPPER_HPP_
