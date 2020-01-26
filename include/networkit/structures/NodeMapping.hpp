/*
 * NodeMapping.hpp
 *
 * Created: 2019-03-28
 * Author: Armin Wiebigke
 */

#ifndef NETWORKIT_STRUCTURES_NODE_MAPPING_HPP_
#define NETWORKIT_STRUCTURES_NODE_MAPPING_HPP_

#include <networkit/graph/Graph.hpp>

namespace NetworKit {

/**
 * Map node IDs between a global graph and a local graph. The mapping is one-to-one.
 */
class NodeMapping {
public:
    NodeMapping() = default;

    explicit NodeMapping(Graph const &G);
    
    explicit NodeMapping(count globalUpperBound);

    /**
     * Add a node to the mapping. Does nothing if the node is already mapped. The local ID of the
     * added node is equal to the current number of added nodes.
     * @param u The node to add.
     * @return true if the node was added, false if node was already mapped
     */
    bool addNode(node u) {
        if (!isMapped(u)) {
            globalToLocal[u] = localToGlobal.size();
            localToGlobal.push_back(u);
            return true;
        }
        return false;
    }

    /**
     * Add a dummy node that has no global node ID. Useful to create a specific mapping from local
     * to global IDs (add dummies for all unused local IDs).
     */
    void addDummy() {
        localToGlobal.push_back(none);
    }

    /**
     * Get the local node from a global node.
     * @param globalNode The global node.
     * @return The mapped local node.
     */
    node toLocal(node globalNode) const {
        assert(globalToLocal[globalNode] != none);
        return globalToLocal[globalNode];
        }

    /**
     * Get the global node from a local node.
     * @param localNode The local node.
     * @return The mapped global node.
     */
    node toGlobal(node localNode) const {
        assert(localToGlobal[localNode] != none);
        return localToGlobal[localNode];
    }

    /**
     * Check if a global node is mapped.
     * @param globalNode The global node.
     * @return True iff the node is mapped.
     */
    bool isMapped(node globalNode) const {
        return globalToLocal[globalNode] != none;
    }

    /**
     * Get the number of mapped nodes.
     * @return The number of mapped nodes.
     */
    count nodeCount() const {
        return localToGlobal.size();
    }

    /**
     * Returns the global node IDs for all mapped nodes.
     * @return A vector of the global nodes.
     */
    const std::vector<node>& globalNodes() const {
        return localToGlobal;
    }

    /**
     * Delete all mappings.
     */
    void reset();

    /**
     * Reset the mapping while keeping only the smallest local nodes.
     * @param end The first local node that will be deleted.
     */
    void reset(index end);

private:
    std::vector<node> globalToLocal;
    std::vector<node> localToGlobal;
};

} /* namespace NetworKit */

#endif // NETWORKIT_STRUCTURES_NODE_MAPPING_HPP_
