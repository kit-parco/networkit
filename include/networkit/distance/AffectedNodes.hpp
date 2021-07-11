// no-networkit-format
/* 
 * File:   AffectedNodes.h
 * Author: paddya
 *
 * Created on 17. Dezember 2016, 16:35
 */

#ifndef NETWORKIT_DISTANCE_AFFECTED_NODES_HPP_
#define NETWORKIT_DISTANCE_AFFECTED_NODES_HPP_

#include <networkit/base/Algorithm.hpp>
#include <networkit/graph/Graph.hpp>
#include <networkit/dynamics/GraphEvent.hpp>

namespace NetworKit {
    
    class AffectedNodes : public Algorithm {
    public:
        /**
         * Constructs the AffectedNodes class for a given graph @a G and a graph event
         * @a event. The run() method computes the set of affected nodes,
         * their distance to the edge modification and, in the case of an edge insertion,
         * an upper bound for the improvement of the harmonic closeness centrality
         * of each affected node.
         * 
         * @param G The graph.
         * @param event The graph event.
         */
        AffectedNodes(const Graph& G, const GraphEvent& event);
        /**
         * Computes the set of affected nodes.
         */
        void run() override;
        /**
         * Returns the set of affected nodes.
         * 
         * @return The set of affected nodes
         */
        std::vector<node> getNodes();
        /**
         * Returns the distances to the edge modification for each node.
         * 
         * @return The distances to the edge modification for each node 
         */
        std::vector<edgeweight> getDistances();
        /**
         * Returns the level-based improvement bounds for each affected node.
         * 
         * @return The improvement upper bound for each affected node, indexed by node ID 
         */
        std::vector<edgeweight> getImprovements();
        
        edgeweight closenessU = 0;
        edgeweight closenessV = 0;
        edgeweight improvementU = 0;
        edgeweight improvementV = 0;
    private:
        /**
         * Handles an edge insertion.
         */
        void addedEdge();
        /**
         * Handles an edge removal.
         */
        void removedEdge();
        /**
         * Runs a complete BFS from @a source while ignoring the edge between 
         * @a source and @a startNeighbor. This is equivalent to running a
         * BFS in G $\\setminus$ (source, startNeighbor).
         * 
         * @param source The source node.
         * @param startNeighbor The neighbor node to be ignored on level 1.
         * @return The distances of each node to @a source.
         */
        std::vector<edgeweight> bfsWithoutStartNeighbor(node source, node startNeighbor);
        
        /**
         * Runs a complete reverse BFS from @a source while ignoring the edge between 
         * @a startNeighbor and @a source. This is equivalent to running a
         * BFS in G - (source, startNeighbor).
         * 
         * @param source The source node.
         * @param startNeighbor The neighbor node to be ignored on level 1.
         * @return The distances of each node to @a source.
         */
        std::vector<edgeweight> reverseBfsWithoutStartNeighbor(node source, node startNeighbor);
        
        /**
         * Runs a BFS from @a source while optionally including an edge between
         * @a source and @a additionalStartNeighbor. This is equivalent to running
         * a BFS from @a source in  G + (source, additionalStartNeighbor).
         * 
         * The algorithm prunes a subtree if it encounters a node @a w for which
         * the new distance is not smaller than the distance in @a distances.
         * 
         * @param source The source node.
         * @param distances The distances without the additional edge.
         * @param additionalStartNeighbor The additional start neighbor.
         * @return The distances of each node to @a source.
         */
        std::pair<std::vector<node>, std::vector<edgeweight>> getAffectedNodes(node source, std::vector<edgeweight>& distances, node additionalStartNeighbor = none);
        
        /**
         * Runs a reverse BFS from @a source while optionally including an edge between
         * @a additionalStartNeighbor and @a source. This is equivalent to running
         * a BFS from @a source in  G + (source, additionalStartNeighbor).
         * 
         * The algorithm prunes a subtree if it encounters a node @a w for which
         * the new distance is not smaller than the distance in @a distances.
         * 
         * @param source The source node.
         * @param distances The distances without the additional edge.
         * @param additionalStartNeighbor The additional start neighbor.
         * @return The distances of each node to @a source.
         */
        std::pair<std::vector<node>, std::vector<edgeweight>> getAffectedNodesBackwards(node source, std::vector<edgeweight>& distances, node additionalStartNeighbor = none);
        const Graph &G;
        const GraphEvent &event;
        std::vector<node> nodes;
        std::vector<edgeweight> distances;
        std::vector<edgeweight> improvements;
    };
}


#endif // NETWORKIT_DISTANCE_AFFECTED_NODES_HPP_

