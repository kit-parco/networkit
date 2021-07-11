// no-networkit-format
/*
 * GMLGraphReader.cpp
 *
 *  Created on: 18.09.2014
 *      Author: Maximilian Vogel (maximilian.vogel@student.kit.edu)
 */

#include <exception>
#include <fstream>
#include <unordered_map>

#include <tlx/define/attribute_fallthrough.hpp>

#include <networkit/auxiliary/Enforce.hpp>
#include <networkit/auxiliary/Log.hpp>
#include <networkit/auxiliary/StringTools.hpp>
#include <networkit/io/GMLGraphReader.hpp>

namespace NetworKit {

Graph GMLGraphReader::read(const std::string& path) {
    std::ifstream graphFile(path);
    Aux::enforceOpened(graphFile);
    std::string line;

    std::unordered_map<std::string,node> nodeMap;

    auto ignoreWhitespaces = [](std::string& s, index i) {
        index newIdx = i;
        index end = s.size();
        while (s[newIdx] == ' ' && newIdx < end) ++newIdx;
        return newIdx;
    };

    auto getPair = [&ignoreWhitespaces](std::string& s) {
        if (s.find("[") < s.size()) {
            throw std::runtime_error("found opening bracket");
        } else if (s.find("]") < s.size()) {
            throw std::runtime_error("found closing bracket");
        }
        index length = s.size();
        index start = ignoreWhitespaces(s,0);
        index i = start;
        while (s[i] != ' ' && i < length) ++i;
        index end = i;
        std::string key = s.substr(start,end-start);
        i = ignoreWhitespaces(s,i+1);
        start = i;
        while (s[i] != ' ' && i < length) ++i;
        end = i;
        std::string value = s.substr(start, end-start);
        return std::make_pair(key,value);
    };

    auto parseNode = [&](Graph& G) {
        bool closed = false;
        if (line.find("node") < line.size() && line.find("[") < line.size()) {
            std::getline(graphFile, line);
            closed = line.find("]") < line.size();
            while (!closed) {
                try {
                    auto pair = getPair(line);
                    if (pair.first == "id") {
                        // currently, node attributes are ignored and only the id is relevant for edges later on
                        node u = G.addNode();
                        nodeMap.insert(std::make_pair(pair.second,u));
                        DEBUG("added node: ",u,", ",pair.second);
                    }
                } catch (const std::exception &e) {
                    return false;
                }
                std::getline(graphFile, line);
                closed = line.find("]") < line.size();
            }
        } else {
            return false;
        }
        if (closed) {
            return true;
        } else {
            return false;
        }
    };

    auto parseEdge = [&](Graph& G) {
        if (line.find("edge") < line.size() && line.find("[") < line.size()) {
            std::getline(graphFile, line);
            node u = 0;
            node v = 0;
            while (!(line.find("]") < line.size())) {
                try {
                    auto pair = getPair(line);
                    if (pair.first == "source") {
                        // currently, node attributes are ignored and only the id is relevant for edges later on
                        u = nodeMap[pair.second];
                    } else if (pair.first == "target") {
                        v = nodeMap[pair.second];
                    }
                } catch (const std::exception &e) {
                    return false;
                }
                std::getline(graphFile,line);
            }
            G.addEdge(u,v);
            DEBUG("added edge ",u ,", ", v);
        } else {
            return false;
        }
        if (line.find("]") < line.size()) {
            return true;
        } else {
            return false;
        }
    };


    auto parseGraph = [&]() {
        Graph G;
        int key_type = 0; // 0 for graph keys, 1 for node, 2 for edges
        bool directed = false;
        std::getline(graphFile, line);
        index pos = line.find("graph");
        bool graphTagFound = pos < line.size();
        if (graphTagFound) {
            if (line.find("[",pos) < line.size()) {
                while (graphFile.good()) {
                    std::getline(graphFile,line);
                    switch (key_type) {
                        case 0:	try {
                                auto pair = getPair(line);
                                // currently, it is only interesting, if the graph is directed
                                if (pair.first == "directed" && pair.second == "1") {
                                    directed = true;
                                    DEBUG("set directed to true");
                                }
                                break;
                            } catch (const std::exception &e) {
                                if (directed) {
                                    G = Graph(0,false,true);
                                } else {
                                    G = Graph(0,false,false);
                                }
                                ++key_type;
                            }
                        case 1: if (parseNode(G)) {
                                break;
                            } else {
                                ++key_type;
                                DEBUG("couldn't parse node; key type is now: ",key_type);
                            }
                        TLX_ATTRIBUTE_FALLTHROUGH;
                        case 2: if (parseEdge(G)) {
                                break;
                            } else {
                                DEBUG("parsing edge went wrong");
                                ++key_type;
                            }
                        TLX_ATTRIBUTE_FALLTHROUGH;
                        default: if (!(line.find("]") < line.size())) {
                                DEBUG("expected closing bracket \"]\", got: ",line);
                             }
                    }
                }// at the end of the file, make sure the closing bracket is there.
            } else {
                throw std::runtime_error("GML graph file broken");
            }
        } else {
            throw std::runtime_error("graph key not found");
        }
        return G;
    };

    Graph G = parseGraph();
    G.shrinkToFit();
    return G;
}

} /* namespace NetworKit */
