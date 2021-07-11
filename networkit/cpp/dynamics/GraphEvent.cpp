// no-networkit-format
/*
 * GraphEvent.cpp
 *
 *  Created on: 02.04.2013
 *      Author: cls
 */

#include <sstream>

#include <networkit/dynamics/GraphEvent.hpp>

namespace NetworKit {

GraphEvent::GraphEvent(GraphEvent::Type type, node u, node v, edgeweight w) : type(type), u(u), v(v), w(w) {}

std::string GraphEvent::toString() const {
    std::stringstream ss;
    if (this->type == GraphEvent::NODE_ADDITION) {
        ss << "an(" << u << ")";
    } else if (this->type == GraphEvent::NODE_REMOVAL) {
        ss << "dn(" << u << ")";
    } else if (this->type == GraphEvent::NODE_RESTORATION) {
        ss << "rn(" << u << ")";
    } else if (this->type == GraphEvent::EDGE_ADDITION) {
        ss << "ae(" << u << "," << v << "," << w << ")";
    } else if (this->type == GraphEvent::EDGE_REMOVAL) {
        ss << "de(" << u << "," << v << ")";
    } else if (this->type == GraphEvent::EDGE_WEIGHT_UPDATE) {
        ss << "ce(" << u << "," << v << ")";
    } else if (this->type == GraphEvent::EDGE_WEIGHT_INCREMENT) {
        ss << "ie(" << u << "," << v << ")";
    } else if (this->type == GraphEvent::TIME_STEP) {
        ss << "st";
    }
    return ss.str();
}

bool GraphEvent::compare(GraphEvent a, GraphEvent b) {
    if (a.type < b.type || (a.type == b.type && a.u < b.u) || (a.type == b.type && a.u == b.u && a.v < b.v) || (a.type == b.type && a.u == b.u && a.v == b.v && a.w < b.w)) return true;
    else return false;
}

bool GraphEvent::equal(GraphEvent a, GraphEvent b) {
    if (a.type == GraphEvent::TIME_STEP && b.type == GraphEvent::TIME_STEP) return true;
    return (a.type == b.type && a.u == b.u && a.v == b.v && a.w == b.w);
}

} /* namespace NetworKit */
