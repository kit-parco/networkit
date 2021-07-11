// no-networkit-format
/*
 * AStar.hpp
 *
 *  Created on: 09.08.2019
 *      Author: Eugenio Angriman <angrimae@hu-berlin.de>
 */

#ifndef NETWORKIT_DISTANCE_A_STAR_HPP_
#define NETWORKIT_DISTANCE_A_STAR_HPP_

#include <networkit/distance/AStarGeneral.hpp>

/**
 * @ingroup distance
 * AStar path-finding algoritm
 */
namespace NetworKit {
class AStar final : public STSP {
public:
    /**
     * Creates the AStarGeneral class for the graph @a G, the source node @a
     * source, and the target node @a target using a template parameter as
     * heuristic function.
     *
     * @param G The graph.
     * @param distanceHeu Vector with a lower bound of the distance from every
     * node to the target.
     * @param source The source node.
     * @param target The target node.
     * @param storePred If true, the algorithm will also store the predecessors
     * and reconstruct a shortest path from @a source and @a target.
     */
    AStar(const Graph &G, const std::vector<double> &distanceHeu, node source, node target,
          bool storePred = true)
        : STSP(G, source, target, storePred),
          astar(AStarGeneral<Heuristic>(G, Heuristic(distanceHeu), source, target, storePred)) {
        if (G.upperNodeIdBound() != distanceHeu.size()) {
            throw std::runtime_error("Error: G.upperNodeIdBound() must be equal to "
                                     "the size of distanceHeu.");
        }
    }

    /*
     * Executes the AStar algorithm.
     */
    void run() override {
        astar.run();
        hasRun = true;
    }

    edgeweight getDistance() const override { return astar.getDistance(); }

    std::vector<node> getPath() const override { return astar.getPath(); }

private:
    struct Heuristic {
    public:
        Heuristic(const std::vector<double> &distanceHeu) : distanceHeu(distanceHeu) {}
        double operator()(node u) const noexcept { return distanceHeu[u]; }

    private:
        const std::vector<double> &distanceHeu;
    };

    AStarGeneral<Heuristic> astar;
};
} // namespace NetworKit

#endif // NETWORKIT_DISTANCE_A_STAR_HPP_
