// no-networkit-format
/*
 * DynamicDorogovtsevMendesGenerator.cpp
 *
 *  Created on: 03.02.2014
 *      Author: cls
 */

#include <networkit/generators/DynamicDorogovtsevMendesGenerator.hpp>

namespace NetworKit {

DynamicDorogovtsevMendesGenerator::DynamicDorogovtsevMendesGenerator() : initial(true), u(0) {}

std::vector<GraphEvent> DynamicDorogovtsevMendesGenerator::generate(count nSteps) {

    std::vector<GraphEvent> stream;

    if (initial) {
        node s1 = u++;
        stream.emplace_back(GraphEvent::NODE_ADDITION, s1);
        node s2 = u++;
        stream.emplace_back(GraphEvent::NODE_ADDITION, s2);
        node s3 = u;
        stream.emplace_back(GraphEvent::NODE_ADDITION, s3);
        edges.emplace_back(s1, s2);
        stream.emplace_back(GraphEvent::EDGE_ADDITION, s1, s2);
        edges.emplace_back(s2, s3);
        stream.emplace_back(GraphEvent::EDGE_ADDITION, s2, s3);
        edges.emplace_back(s3, s1);
        stream.emplace_back(GraphEvent::EDGE_ADDITION, s3, s1);

        stream.emplace_back(GraphEvent::TIME_STEP);
        initial = false;
    }

    for (index i = 0; i < nSteps; ++i) {
        ++u; // new node
        stream.emplace_back(GraphEvent::NODE_ADDITION, u);
        // select random edge
        index e = Aux::Random::integer(edges.size() - 1);
        node s = edges[e].first;
        node t = edges[e].second;
        edges.emplace_back(s, u);
        edges.emplace_back(t, u);
        // connect node
        stream.emplace_back(GraphEvent::EDGE_ADDITION, u, s);
        stream.emplace_back(GraphEvent::EDGE_ADDITION, u, t);

        stream.emplace_back(GraphEvent::TIME_STEP);
    }

    return stream;
}

} /* namespace NetworKit */
