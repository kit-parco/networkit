/*
 * Curveball.cpp
 *
 *  Created on: 26.05.2018
 *      Author:  Hung Tran <htran@ae.cs.uni-frankfurt.de>, Manuel Penschuck <networkit@manuel.jetzt>
 */

#include "CurveballImpl.hpp"
#include <networkit/auxiliary/Random.hpp>
#include <networkit/randomization/Curveball.hpp>

namespace NetworKit {

Curveball::Curveball(const Graph &G) : impl(new CurveballDetails::CurveballIM{G}) {}

// We have to define a "default" destructor here, since the definition of
// CurveballDetails::CurveballImpl is not known in the header file
Curveball::~Curveball() = default;

void Curveball::run(const CurveballDetails::trade_vector &trades) {
    impl->run(trades);
}

Graph Curveball::getGraph(bool parallel) {
    return impl->getGraph(parallel);
}

std::string Curveball::toString() const {
    return "Curveball";
}

count Curveball::getNumberOfAffectedEdges() const {
    return impl->getNumberOfAffectedEdges();
}

} // namespace NetworKit
