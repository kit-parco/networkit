// no-networkit-format
/*
 * DynamicGenerator.cpp
 *
 *  Created on: 03.04.2013
 *      Author: cls
 */

#include <networkit/generators/DynamicGraphSource.hpp>

namespace NetworKit {

DynamicGraphSource::DynamicGraphSource() : Gproxy(NULL), G(NULL), graphSet(false), graphInitialized(false) {
    // Graph and GraphEventProxy are set by calling newGraph
}

void DynamicGraphSource::generateWhile(std::function<bool(void)> cont) {
    while (cont()) {
        this->generate();
    }
}

void DynamicGraphSource::generateNodes(count n) {
    while (G->numberOfNodes() < n) {
        this->generate();
    }
}

void DynamicGraphSource::generateEdges(count m) {
    while (G->numberOfEdges() < m) {
        this->generate();
    }
}

GraphEventProxy* DynamicGraphSource::newGraph() {
    this->G = new Graph(0);
    this->Gproxy = new GraphEventProxy(*(this->G));
    // not returning proxy because only generator needs write access to graph
    this->graphSet = true;
    return this->Gproxy;
}

void DynamicGraphSource::generateTimeSteps(count t) {
    while (G->time() < t) {
        this->generate();
    }
    G->shrinkToFit(); // TODO shrinkToFit: is this method supposed be be call often?
}

} /* namespace NetworKit */
