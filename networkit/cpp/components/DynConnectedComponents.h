/*
* DynConnectedComponents.cpp
*
*  Created on: June 2017
*      Author: Eugenio Angriman
*/

#ifndef DYNCONNECTEDCOMPONENTS_H_
#define DYNCONNECTEDCOMPONENTS_H_

#include "../graph/Graph.h"
#include "../base/Algorithm.h"
#include "../base/DynAlgorithm.h"
#include "../dynamics/GraphEvent.h"

namespace NetworKit {

	/**
	* @ingroup components
	* Determines and updates the connected components of an undirected graph.
	*/
	class DynConnectedComponents : public Algorithm, public DynAlgorithm {

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
		const Graph& G;
		std::vector<bool> isTree;
		std::vector<index> components;
		std::map<index, count> compSize;
		std::map<std::pair<node, node>, index> edgesMap;
		std::vector<count> tmpDistances;
		std::queue<index> componentIds;
		bool distancesInit;
		bool hasRun;
	};

	inline count DynConnectedComponents::componentOfNode(node u) {
		assert(u <= G.upperNodeIdBound());
		assert (components[u] != none);
		if (!hasRun) throw std::runtime_error("run method has not been called");
		return components[u];
	}

	inline count DynConnectedComponents::numberOfComponents() {
		if (!hasRun) throw std::runtime_error("run method has not been called");
		return this->compSize.size();
	}

	inline std::map<index, count> DynConnectedComponents::getComponentSizes() {
		if (!hasRun){
			throw std::runtime_error("run method has not been called");
		}
		return compSize;
	}

}

#endif /* DYNCONNECTEDCOMPONENTS_H_ */
