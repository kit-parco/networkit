// no-networkit-format
/*
 * LayoutAlgorithm.hpp
 *
 *  Created on: May 20 2015
 *      Author: Christian Staudt
 */

#ifndef NETWORKIT_LAYOUT_LAYOUT_ALGORITHM_HPP_
#define NETWORKIT_LAYOUT_LAYOUT_ALGORITHM_HPP_

#include <networkit/graph/Graph.hpp>

namespace NetworKit {

/**
 * @ingroup layout
 *
 * Base class for graph layout algorithms, i.e. algorithms that compute
 * a layout in form of 2D coordinates for nodes.
 */
class LayoutAlgorithm {


public:

    LayoutAlgorithm(const Graph& G) : G(G) {

    };

    virtual void run() = 0;

    virtual std::vector<std::pair<double, double>> getLayout() {
        return layout;
    };

protected:

    const Graph* G;
    std::vector<std::pair<double, double>> layout;

};

} /* namespace NetworKit */
#endif // NETWORKIT_LAYOUT_LAYOUT_ALGORITHM_HPP_
