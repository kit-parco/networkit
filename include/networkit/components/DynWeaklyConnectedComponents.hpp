// no-networkit-format
/*
* DynWeaklyConnectedComponents.h
*
*  Created on: June 20, 2017
*      Author: Eugenio Angriman
*/

#ifndef NETWORKIT_COMPONENTS_DYN_WEAKLY_CONNECTED_COMPONENTS_HPP_
#define NETWORKIT_COMPONENTS_DYN_WEAKLY_CONNECTED_COMPONENTS_HPP_

#include <map>
#include <vector>

#include <networkit/graph/Graph.hpp>
#include <networkit/base/Algorithm.hpp>
#include <networkit/base/DynAlgorithm.hpp>
#include <networkit/dynamics/GraphEvent.hpp>

namespace NetworKit {

    /**
    * @ingroup components
    * Determines and updates the weakly connected components of a directed
    * graph.
    */
    class DynWeaklyConnectedComponents final : public Algorithm, public DynAlgorithm {

    public:
        /**
        * Create DynWeaklyConnectedComponents class for Graph @a G.
        *
        * @param G The graph.
        */
        DynWeaklyConnectedComponents(const Graph& G);

        /**
        * This method determines the weakly connected components for the graph
        * given in the constructor.
        */
        void run() override;

        /**
        * Updates the weakly connected components after an edge insertion or
        * deletion.
        *
        * @param[in]	event	The event that happened (edge insertion
        *                       or deletion).
        */
        void update(GraphEvent event) override;

        /**
        * Updates the weakly connected components after a batch of edge events
        * (insertions or deletions).
        *
        * @param[in] batch  A vector that contains a batch of edge events
        *                   (insertions or deletions).
        */
        void updateBatch(const std::vector<GraphEvent>& batch) override;

        /**
        * Get the number of weakly connected components.
        *
        * @return The number of weakly connected components.
        */
        count numberOfComponents();

        /**
        * Get the the component in which node @a u is.
        *
        * @param[in]	u	The node.
        */
        count componentOfNode(node u);

        /**
        * @return Map from component to size.
        */
        std::map<index, count> getComponentSizes();

        /**
        * @return Vector of components, each stored as (unordered) set of nodes.
        */
        std::vector<std::vector<node>> getComponents();


    private:
        void init();
        void addEdge(node u, node v);
        void removeEdge(node u, node v);
        void updateComponent(index c, node w, std::queue<node>& q, node v);
        void indexEdges();
        void insertEdgeIntoMap(node u, node v, edgeid eid);
        edgeid updateMapAfterAddition(node u, node v);
        void updateTreeAfterAddition(edgeid eid, bool partOfTree);
        void reverseBFS(node u, node v);
        index getEdgeId(node u, node v);
        index nextAvailableComponentId(bool eraseId = true);
        bool visitNodeReversed(
            node u,
            node s,
            node w,
            node v,
            count d,
            std::queue<node>& q,
            bool& nextEdgeFound,
            count level
        );
        std::pair<node, node> makePair(node u, node v);

        // Pointer to the graph
        const Graph* G;

        // This vector associates each node to a component
        std::vector<index> components;

        // Whether each edge is part of a spanning tree. If an edge is removed
        // then a components could split.
        std::vector<bool> isTree;

        // Maps the edges to their IDs
        std::map<std::pair<node, node>, edgeid> edgesMap;

        // Temporary distances used for BFSs after edge updates.
        std::vector<count> tmpDistances;

        // Stores all components <Key: component ID, Value: component size>
        std::map<index, count> compSize;

        // Stores available component IDs.
        std::queue<index> componentIds;
    };


    inline count DynWeaklyConnectedComponents::componentOfNode(node u) {
        assert (components[u] != none);
        assureFinished();
        return components[u];
    }

    inline count DynWeaklyConnectedComponents::numberOfComponents() {
        assureFinished();
        return this->compSize.size();
    }

    inline std::map<index, count> DynWeaklyConnectedComponents::getComponentSizes() {
        assureFinished();
        return compSize;
    }

}


#endif // NETWORKIT_COMPONENTS_DYN_WEAKLY_CONNECTED_COMPONENTS_HPP_
