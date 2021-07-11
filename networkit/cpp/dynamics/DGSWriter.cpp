// no-networkit-format
/*
 * DGSWriter.cpp
 *
 *  Created on: 14.01.2014
 *      Author: cls
 */

#include <fstream>

#include <networkit/dynamics/DGSWriter.hpp>

namespace NetworKit {

void DGSWriter::write(std::vector<GraphEvent>& stream, const std::string& path) {
    std::ofstream out(path);

    out << "DGS004\nnoname 0 0\n";

    for (const auto &ev : stream) {
        switch (ev.type) {
            case GraphEvent::NODE_ADDITION : {
                out << "an " << ev.u << "\n";
                break;
            }
            case GraphEvent::NODE_REMOVAL : {
                out << "dn " << ev.u << "\n";
                break;
            }
            case GraphEvent::NODE_RESTORATION : {
                out << "rn" << ev.u << "\n";
                break;
            }
            case GraphEvent::EDGE_ADDITION : {
                out << "ae " << ev.u << "-" << ev.v << " " << ev.u << " " << ev.v << " weight=" << ev.w << "\n";
                break;
            }
            case GraphEvent::EDGE_REMOVAL : {
                out << "de " << ev.u << "-" << ev.v << "\n";
                break;
            }
            case GraphEvent::EDGE_WEIGHT_UPDATE : {
                out << "ce " << ev.u << "-" << ev.v  << " weight=" << ev.w << "\n";
                break;
            }
            case GraphEvent::EDGE_WEIGHT_INCREMENT : {
                out << "ie " << ev.u << "-" << ev.v  << " weight=" << ev.w << "\n";
                break;
            }
            case GraphEvent::TIME_STEP : {
                out << "st" << "\n";
                break;
            }
            default: {
                throw std::runtime_error("unknown event type");
            }
        }
    }
}

} /* namespace NetworKit */
