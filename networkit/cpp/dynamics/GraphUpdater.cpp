// no-networkit-format
/*
 * GraphUpdater.cpp
 *
 *  Created on: 27.12.2013
 *      Author: cls
 */

#include <networkit/auxiliary/Log.hpp>
#include <networkit/dynamics/GraphUpdater.hpp>

namespace NetworKit {

GraphUpdater::GraphUpdater(Graph& G) : G(&G) {}

void GraphUpdater::update(const std::vector<GraphEvent>& stream) {
    for (GraphEvent ev : stream) {
        TRACE("event: " , ev.toString());
        switch (ev.type) {
            case GraphEvent::NODE_ADDITION : {
                G->addNode();
                break;
            }
            case GraphEvent::NODE_REMOVAL : {
                G->removeNode(ev.u);
                break;
            }
            case GraphEvent::NODE_RESTORATION :{
                G->restoreNode(ev.u);
                break;
            }
            case GraphEvent::EDGE_ADDITION : {
                G->addEdge(ev.u, ev.v, ev.w);
                break;
            }
            case GraphEvent::EDGE_REMOVAL : {
                G->removeEdge(ev.u, ev.v);
                break;
            }
            case GraphEvent::EDGE_WEIGHT_UPDATE : {
                G->setWeight(ev.u, ev.v, ev.w);
                break;
            }
            case GraphEvent::EDGE_WEIGHT_INCREMENT : {
                G->setWeight(ev.u, ev.v, G->weight(ev.u, ev.v) + ev.w);
                break;
            }
            case GraphEvent::TIME_STEP : {
                G->timeStep();
                break;
            }
            default: {
                throw std::runtime_error("unknown event type");
            }
        }
    }
    // record graph size
    size.emplace_back(G->numberOfNodes(), G->numberOfEdges());
}

std::vector<std::pair<count, count>> GraphUpdater::getSizeTimeline() {
    return size;
}

} /* namespace NetworKit */
