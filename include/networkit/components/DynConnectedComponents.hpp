// no-networkit-format
/*
* DynConnectedComponents.hpp
*
*  Created on: June 2017
*      Author: Eugenio Angriman
*/

#ifndef NETWORKIT_COMPONENTS_DYN_CONNECTED_COMPONENTS_HPP_
#define NETWORKIT_COMPONENTS_DYN_CONNECTED_COMPONENTS_HPP_

#include <map>
#include <vector>
#include <queue>

#include <networkit/graph/Graph.hpp>
#include <networkit/base/Algorithm.hpp>
#include <networkit/base/DynAlgorithm.hpp>
#include <networkit/dynamics/GraphEvent.hpp>

namespace NetworKit {

    /**
    * @ingroup components
    * Determines and updates the connected components of an undirected graph.
    */
    class DynConnectedComponents final : public Algorithm, public DynAlgorithm {

    public:
        /**
        * Create ConnectedComponents class for Graph @a G.
        *
        * @param G The graph.
        */
        DynConnectedComponents(const Graph& G);

        /**
        * This method determines the connected components for the graph given in
        *  the constructor.
        */
        void run() override;

        /**
        * Updates the connected components after an edge insertion or deletion.
        *
        * @param[in]	event	The event that happened (edge deletion or
        * insertion).
        */
        void update(GraphEvent e) override;

        /**
        * Updates the connected components after a batch of edge insertions or
        * deletions.
        *
        * @param[in] batch	A vector that contains a batch of edge insertions or
        *					deletions.
        */
        void updateBatch(const std::vector<GraphEvent>& batch) override;

        /**
        * Get the number of connected components.
        *
        * @return The number of connected components.
        */
        count numberOfComponents();

        /**
        * Returns the the component in which node @a u is.
        *
        * @param[in]	u	The node.
        */
        count componentOfNode(node u);

        /**
        * Returns the map from component to size.
        */
        std::map<index, count> getComponentSizes();

        /**
        * @return Vector of components, each stored as (unordered) set of nodes.
        */
        std::vector<std::vector<node>> getComponents();

    private:
        void addEdge(node u, node v);
        void removeEdge(node u, node v);
        void addEdgeDirected(node u, node v);
        void removeEdgeDirected(node u, node v);
        void reverseBFS(node u, node v);
        index nextAvailableComponentId(bool eraseId = true);
        void indexEdges();
        void insertEdgeIntoMap(node u, node v, edgeid eid);
        index getEdgeId(node u, node v);
        // Returns true and the corresponding edge id if the new edge was not
        // into the original graph.
        std::pair<bool, edgeid> updateMapAfterAddition(node u, node v);
        void init();
        std::pair<node, node> makePair(node u, node v);
        const Graph* G;
        std::vector<bool> isTree;
        std::vector<index> components;
        std::map<index, count> compSize;
        std::map<std::pair<node, node>, index> edgesMap;
        std::vector<count> tmpDistances;
        std::queue<index> componentIds;
        bool distancesInit;
    };

    inline count DynConnectedComponents::componentOfNode(node u) {
        assert(u <= G->upperNodeIdBound());
        assert (components[u] != none);
        assureFinished();
        return components[u];
    }

    inline count DynConnectedComponents::numberOfComponents() {
        assureFinished();
        return this->compSize.size();
    }

    inline std::map<index, count> DynConnectedComponents::getComponentSizes() {
        assureFinished();
        return compSize;
    }

}

#endif // NETWORKIT_COMPONENTS_DYN_CONNECTED_COMPONENTS_HPP_
