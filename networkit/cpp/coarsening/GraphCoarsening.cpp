/*
 * GraphCoarsening.cpp
 *
 *  Created on: 30.10.2012
 *      Author: Christian Staudt (christian.staudt@kit.edu)
 */

#include <networkit/coarsening/GraphCoarsening.hpp>

namespace NetworKit {

GraphCoarsening::GraphCoarsening(const Graph& G) : Algorithm(), G(&G) {

}

const Graph& GraphCoarsening::getCoarseGraph() const {
    if(!hasRun) {
        throw std::runtime_error("Call run()-method first.");
    }
    return Gcoarsened;
}

const std::vector<node>& GraphCoarsening::getFineToCoarseNodeMapping() const {
    if(!hasRun) {
        throw std::runtime_error("Call run()-method first.");
    }
    return nodeMapping;
}


std::map<node, std::vector<node> > GraphCoarsening::getCoarseToFineNodeMapping() const {
    assureFinished();

    std::map<node, std::vector<node>> reverseMap;
    Gcoarsened.forNodes([&](node v_){
        std::vector<node> empty;
        reverseMap[v_] = empty;
    });


    G->forNodes([&](node v) {
        node v_ = nodeMapping[v];
        reverseMap[v_].push_back(v);
    });

    return reverseMap;
}


}
