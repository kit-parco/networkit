// no-networkit-format
/*
 * GraphUpdater.hpp
 *
 *  Created on: 27.12.2013
 *      Author: cls
 */

#ifndef NETWORKIT_DYNAMICS_GRAPH_UPDATER_HPP_
#define NETWORKIT_DYNAMICS_GRAPH_UPDATER_HPP_

#include <networkit/dynamics/GraphEvent.hpp>
#include <networkit/graph/Graph.hpp>

namespace NetworKit {

/**
 * @ingroup dynamics
 */
class GraphUpdater final {

public:

    GraphUpdater(Graph& G);

    void update(const std::vector<GraphEvent>& stream);

    std::vector<std::pair<count, count> > getSizeTimeline();

    static bool compare(GraphEvent a, GraphEvent b);
    static bool equal(GraphEvent a, GraphEvent b);

private:

    Graph* G;
    std::vector<std::pair<count, count>> size;
};

} /* namespace NetworKit */

#endif // NETWORKIT_DYNAMICS_GRAPH_UPDATER_HPP_
