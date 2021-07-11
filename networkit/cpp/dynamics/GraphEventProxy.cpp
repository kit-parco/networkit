// no-networkit-format
/*
 * GraphEventProxy.cpp
 *
 *  Created on: 03.04.2013
 *      Author: cls
 */

#include <networkit/auxiliary/Log.hpp>
#include <networkit/dynamics/GraphEventProxy.hpp>

namespace NetworKit {

GraphEventProxy::GraphEventProxy() {}

GraphEventProxy::GraphEventProxy(Graph& G) {
    this->G = &G;
}

node GraphEventProxy::addNode() {
    node u = this->G->addNode();
    for (GraphEventHandler* observer : this->observers) {
        observer->onNodeAddition(u);
    }
    return u;
}

void GraphEventProxy::removeNode(node u) {
    this->G->removeNode(u);
    for (GraphEventHandler* observer : this->observers) {
        observer->onNodeRemoval(u);
    }
}

void GraphEventProxy::restoreNode(node u) {
    this->G->restoreNode(u);
    for (GraphEventHandler* observer : this->observers) {
        observer->onNodeRestoration(u);
    }
}
void GraphEventProxy::addEdge(node u, node v, edgeweight weight) {
    this->G->addEdge(u, v, weight);
    for (GraphEventHandler* observer : this->observers) {
        observer->onEdgeAddition(u, v);
    }
}

void GraphEventProxy::removeEdge(node u, node v) {
    this->G->removeEdge(u, v);
    for (GraphEventHandler* observer : this->observers) {
        observer->onEdgeRemoval(u, v);
    }
}


void GraphEventProxy::setWeight(node u, node v, edgeweight w) {
    edgeweight wOld = this->G->weight(u, v);
    this->G->setWeight(u, v, w);
    for (GraphEventHandler* observer : this->observers) {
        observer->onWeightUpdate(u, v, wOld, w);
    }
}

void GraphEventProxy::incrementWeight(node u, node v, edgeweight delta) {
    edgeweight wOld = this->G->weight(u, v);
    this->G->setWeight(u, v, wOld+delta);
    for (GraphEventHandler* observer : this->observers) {
        observer->onWeightIncrement(u, v, wOld, delta);
    }
}

void GraphEventProxy::timeStep() {
    WARN("GraphEventProxy::timeStep is deprecated and will not be supported in future releases.");
    // increment time step counter in G
    this->G->timeStep();
    for (GraphEventHandler* observer : this->observers) {
        observer->onTimeStep();
    }
}

void GraphEventProxy::registerObserver(GraphEventHandler* observer) {
    this->observers.push_back(observer);
}

} /* namespace NetworKit */
