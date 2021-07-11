// no-networkit-format
/*
 * METISGraphWriter.cpp
 *
 *  Created on: 30.01.2013
 *      Author: Christian Staudt
 */

#include <fstream>

#include <networkit/auxiliary/Enforce.hpp>
#include <networkit/graph/GraphTools.hpp>
#include <networkit/io/METISGraphWriter.hpp>

namespace NetworKit {

void METISGraphWriter::write(const Graph &G, const std::string &path) {
    this->write(G, G.isWeighted(), path);
}

void METISGraphWriter::write(const Graph &G, bool weighted, const std::string &path) {
    if (G.isDirected()) {
        throw std::invalid_argument{"METIS does not support directed graphs"};
    }
    std::ofstream file(path);
    Aux::enforceOpened(file);

    count n = G.numberOfNodes();
    count m = G.numberOfEdges();

    file << n << " " << m << " " << int{weighted} << '\n';
    if (G.numberOfNodes() != G.upperNodeIdBound()) {
        auto nodeIds = GraphTools::getContinuousNodeIds(G);
        if (weighted == false) {
            G.forNodes([&](node u) {
                G.forNeighborsOf(u, [&](node v){
                    file << nodeIds[v] + 1 << " ";
                });
                file << '\n';
            });
        } else {
                G.forNodes([&](node u) {
                G.forNeighborsOf(u, [&](node v, edgeweight w){
                    file << nodeIds[v] + 1 << " " << w <<"\t";
                });
                file << '\n';
            });
        }
    } else {
        if (weighted == false) {
            G.forNodes([&](node u) {
                G.forNeighborsOf(u, [&](node v){
                    file << v + 1 << " ";
                });
                file << '\n';
            });
        } else {
                G.forNodes([&](node u) {
                G.forNeighborsOf(u, [&](node v, edgeweight w){
                    file << v + 1 << " " << w <<"\t";
                });
                file << '\n';
            });
        }
    }

}

} /* namespace NetworKit */
