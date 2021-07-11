// no-networkit-format
/*
 * Sfigality.h
 *
 *  Created on: 20.01.2016
 *      Author:Christian Staudt
 */

#ifndef NETWORKIT_CENTRALITY_SFIGALITY_HPP_
#define NETWORKIT_CENTRALITY_SFIGALITY_HPP_

#include <networkit/centrality/Centrality.hpp>

namespace NetworKit {

/**
 * @ingroup centrality
 * The sfigality of a node is the ratio of neighboring nodes that have
 * a higher degree than the node itself.
 */
class Sfigality: public Centrality {
public:
    /**
     * Constructs the Sfigality class for the given Graph @a G. Sfigality is a new type of
     * node centrality measures that is high if neighboring nodes have a higher degree, e.g. in social networks, if your friends have more friends than you. Formally:
     *
     * $$\\sigma(u) = \\frac{| \\{ v: \\{u,v\\} \\in E, deg(u) < deg(v) \\} |}{ deg(u) }$$
     *
     * @param G The graph.

     */
    Sfigality(const Graph& G);

    void run() override;

    /**
     * Not implemented.
     * The maximum sfigality is 1, when all neighbors of a node have a higher degree.
     */
    double maximum() override;
};

} /* namespace NetworKit */

#endif // NETWORKIT_CENTRALITY_SFIGALITY_HPP_
