/*
 * ClusteringFunctionFactory.hpp
 *
 * Created on: 2019-11-12
 * Author: Armin Wiebigke
  */

#ifndef NETWORKIT_COMMUNITY_CLUSTERING_FUNCTION_FACTORY_HPP_
#define NETWORKIT_COMMUNITY_CLUSTERING_FUNCTION_FACTORY_HPP_

#include <networkit/graph/Graph.hpp>
#include <networkit/structures/Partition.hpp>

namespace NetworKit {

using ClusteringFunction = std::function<Partition(const Graph &)>;

class ClusteringFunctionFactory {
public:
    virtual ClusteringFunction getFunction() const = 0;
};

} // namespace NetworKit

#endif // NETWORKIT_COMMUNITY_CLUSTERING_FUNCTION_FACTORY_HPP_
