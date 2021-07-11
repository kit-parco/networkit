// no-networkit-format
/*
 * METISGraphReader.cpp
 *
 *  Created on: 17.01.2013
 *      Author: Christian Staudt
 */

#include <networkit/auxiliary/Enforce.hpp>
#include <networkit/auxiliary/Log.hpp>
#include <networkit/auxiliary/StringTools.hpp>
#include <networkit/graph/GraphBuilder.hpp>
#include <networkit/io/METISGraphReader.hpp>
#include <networkit/io/METISParser.hpp>

namespace NetworKit {

Graph METISGraphReader::read(const std::string& path) {

    METISParser parser(path);

    std::tuple<count, count, index, count> header = parser.getHeader();
    count n = std::get<0>(header);
    count m = std::get<1>(header);
    index fmt = std::get<2>(header);
    count ncon = std::get<3>(header);

    bool weighted;
    if (fmt % 10 == 1) {
        weighted = true;
        DEBUG("graph has been identified as weighted");
    } else {
        weighted = false;
    }
    count ignoreFirst = 0;
    if (fmt / 10 == 1) {
        DEBUG("first ",ncon," value(s) will be ignored");
        ignoreFirst = ncon;
    }

    GraphBuilder b(n, weighted);
    std::string graphName = Aux::StringTools::split(Aux::StringTools::split(path, '/').back(), '.').front();

    INFO("\n[BEGIN] reading graph G(n=", n, ", m=", m, ") from METIS file: ", graphName);

#ifndef NETWORKIT_RELEASE_LOGGING
    double p = 0.0; // percentage for progress bar
#endif
    node u = 0; // begin with 0
    count edgeCounter = 0;
    if (weighted == 0) {
        while (parser.hasNext() && u < n) {
            std::vector<node> adjacencies = parser.getNext(ignoreFirst);
            edgeCounter += adjacencies.size();
            for (index i=0; i < adjacencies.size(); i++) {
                if (adjacencies[i] == 0) {
                    ERROR("METIS Node ID should not be 0, edge ignored.");
                    continue;
                }
                Aux::Checkers::Enforcer::enforce(adjacencies[i] > 0 && adjacencies[i] <= n);
                node v = adjacencies[i] - 1;// METIS-indices are 1-based
                // correct edgeCounter for selfloops
                edgeCounter += (u == v);
                b.addHalfEdge(u, v);
            }
            u++; // next node
#ifndef NETWORKIT_RELEASE_LOGGING
            if ((u % ((n + 10)/10)) == 0) {
                p = ((double) (u-1) / (double) n) * 100;
                DEBUG(p, "% ");
            }
#endif
        }
    } else {
        while (parser.hasNext() && u < n) {
            std::vector<std::pair<node,double>> adjacencies = parser.getNextWithWeights(ignoreFirst);
            edgeCounter += adjacencies.size();
            DEBUG("node ",u," has ",adjacencies.size()," edges");
            for (index i=0; i < adjacencies.size(); i++) {
                if (adjacencies[i].first == 0) {
                    ERROR("METIS Node ID should not be 0, edge ignored.");
                    continue;
                }
                Aux::Checkers::Enforcer::enforce(adjacencies[i].first > 0 && adjacencies[i].first <= n);
                node v = adjacencies[i].first - 1; 	// METIS-indices are 1-based
                // correct edgeCounter for selfloops
                edgeCounter += (u == v);
                double weight = adjacencies[i].second;
                b.addHalfEdge(u, v, weight);
                TRACE("(",u,",",v,",",adjacencies[i].second,")");
            }
            u += 1; // next node
#ifndef NETWORKIT_RELEASE_LOGGING
            if ((u % ((n + 10)/10)) == 0) {
                p = ((double) (u-1) / (double) n) * 100;
                DEBUG(p, "% ");
            }
#endif
        }
    }

    auto G = b.toGraph(false);

    if (G.numberOfEdges() != m) {
        ERROR("METIS file ", path," is corrupted: actual number of added edges doesn't match the specifed number of edges");
    }
    if (edgeCounter != 2 * m) {
        WARN("METIS file is corrupted: not every edge is listed twice");
    }

    INFO("\n[DONE]\n");
    return G;
}

} /* namespace NetworKit */
