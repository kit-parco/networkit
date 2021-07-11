// no-networkit-format
/*
 * GlobalClusteringCoefficient.h
 *
 *  Created on: 12.11.2013
 */

#ifndef NETWORKIT_GLOBAL_GLOBAL_CLUSTERING_COEFFICIENT_HPP_
#define NETWORKIT_GLOBAL_GLOBAL_CLUSTERING_COEFFICIENT_HPP_

#include <networkit/graph/Graph.hpp>

namespace NetworKit {

/**
 * @ingroup global
 */
class GlobalClusteringCoefficient {

public:  
    virtual double approximate(const Graph& G, int k);
};

} /* namespace NetworKit */
#endif // NETWORKIT_GLOBAL_GLOBAL_CLUSTERING_COEFFICIENT_HPP_
