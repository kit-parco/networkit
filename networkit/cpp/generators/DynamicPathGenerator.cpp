// no-networkit-format
/*
 * DynamicPathGenerator.cpp
 *
 *  Created on: 14.01.2014
 *      Author: cls
 */

#include <networkit/generators/DynamicPathGenerator.hpp>

namespace NetworKit {

std::vector<GraphEvent> DynamicPathGenerator::generate(count nSteps) {

    std::vector<GraphEvent> stream;
    node u = G.addNode();
    stream.emplace_back(GraphEvent::NODE_ADDITION, u);

    count step = 0;
    while (step < nSteps) {
        node v = G.addNode();
        stream.emplace_back(GraphEvent::NODE_ADDITION, v);
        G.addEdge(u, v);
        stream.emplace_back(GraphEvent::EDGE_ADDITION, u, v, 1.0);
        u = v;
        stream.emplace_back(GraphEvent::TIME_STEP);
        step += 1;
    }
    return stream;
}

} /* namespace NetworKit */
