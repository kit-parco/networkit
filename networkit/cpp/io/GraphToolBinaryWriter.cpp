// no-networkit-format
/*
 * GraphToolBinaryReader.cpp
 *
 *  Created on: 02.12.14
 *      Author: Maximilian Vogel
 */

#include <fstream>

#include <networkit/auxiliary/Log.hpp>
#include <networkit/auxiliary/Enforce.hpp>
#include <networkit/graph/GraphTools.hpp>
#include <networkit/io/GraphToolBinaryWriter.hpp>

namespace NetworKit {
GraphToolBinaryWriter::GraphToolBinaryWriter(bool littleEndianness) : littleEndianness(littleEndianness) {}

void GraphToolBinaryWriter::write(const Graph &G, const std::string &path) {
    std::ofstream file(path, std::ios::binary | std::ios::out);
    Aux::enforceOpened(file);
    writeHeader(file);
    writeComment(file);
    // write directed byte
    char directed = (char)G.isDirected();
    file.write(&directed,1);
    // write number of nodes to the file
    uint64_t n = G.numberOfNodes();
    writeType<uint64_t>(file,8,n);
    // write the adjacencies to the file
    writeAdjacencies(file,G);
    // close the file
    file.close();
}

uint8_t GraphToolBinaryWriter::getAdjacencyWidth(uint64_t n) {
    if (n < (uint64_t)1 << 8) {
        return 1;
    } else if (n < (uint64_t)1 << 16) {
        return 2;
    } else if (n < (uint64_t)1 << 32) {
        return 4;
    } else {
        return 8;
    } // error handling?
}

void GraphToolBinaryWriter::writeHeader(std::ofstream &file) {
    uint8_t header[8] = {0xe2, 0x9b, 0xbe, 0x20, 0x67, 0x74, 0x01, 0x00};
    header[7] |= (uint8_t) !this->littleEndianness;
    file.write((char*)header,8);
}

void GraphToolBinaryWriter::writeComment(std::ofstream &file) {
    std::string s = "";
    uint64_t size = (uint64_t)s.size();
    writeType<uint64_t>(file,8,size);
    if (size > 0) {
        file.write(s.c_str(),size);
    }
}

void GraphToolBinaryWriter::writeAdjacencies(std::ofstream &file, const Graph &G) {
    // value of numNodes determines the size of the unsigned integer type storing the node ids
    int width = (int)getAdjacencyWidth(G.numberOfNodes());
    // if the node ids aren't continuous, a node id map becomes necessary
    if (G.upperNodeIdBound() == G.numberOfNodes()) {
        if (!G.isDirected()) {
            // iterate over each node
            G.forNodes([&](node u){
                std::vector<node> adjacencies;
                // collect number of adjacencies for undirected case
                G.forNeighborsOf(u,[&](node v){
                    if (v <= u)
                        adjacencies.push_back(v);
                });
                // write number of adjacencies to the file
                count deg = adjacencies.size();
                writeType<uint64_t>(file,8,deg);
                // write each adjacent to the file
                for(auto v : adjacencies) {
                    writeType<uint64_t>(file,width,v);
                }
            });
        } else {
            // iterate over each node
            G.forNodes([&](node u){
                // write number of adjacencies to the file
                count deg = G.degree(u);
                writeType<uint64_t>(file,8,deg);
                // write each adjacent to the file
                G.forNeighborsOf(u,[&](node v){
                    writeType<uint64_t>(file,width,v);
                });
            });
        }
    } else {
        auto nodeIdMap = GraphTools::getContinuousNodeIds(G);
        if (!G.isDirected()) {
            // iterate over each node
            G.forNodes([&](node u){
                std::vector<node> adjacencies;
                // collect number of adjacencies for undirected case
                G.forNeighborsOf(u,[&](node v){
                    if (v <= u)
                        adjacencies.push_back(v);
                });
                // write number of adjacencies to the file
                count deg = adjacencies.size();
                writeType<uint64_t>(file,8,deg);
                // write each adjacent to the file
                for(auto& v : adjacencies) {
                    writeType<uint64_t>(file,width,nodeIdMap[v]);
                }
            });
        } else {
            // iterate over each node
            G.forNodes([&](node u){
                // write number of adjacencies to the file
                count deg = G.degree(u);
                writeType<uint64_t>(file,8,deg);
                // write each adjacent to the file
                G.forNeighborsOf(u,[&](node v){
                    writeType<uint64_t>(file,width,nodeIdMap[v]);
                });
            });
        }
    }
}
} // namespace NetworKit
