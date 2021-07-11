/*
 * EdgeListReader.cpp
 *
 *  Created on: 18.06.2013
 *      Author: cls
 */

#include <fstream>
#include <sstream>

#include <networkit/auxiliary/Enforce.hpp>
#include <networkit/auxiliary/Log.hpp>
#include <networkit/graph/GraphTools.hpp>
#include <networkit/io/EdgeListReader.hpp>
#include <networkit/io/MemoryMappedFile.hpp>

namespace NetworKit {

EdgeListReader::EdgeListReader(char separator, node firstNode, const std::string &commentPrefix,
                               bool continuous, bool directed)
    : separator(separator), commentPrefix(commentPrefix), firstNode(firstNode),
      continuous(continuous), mapNodeIds(), directed(directed) {
    if (!continuous && firstNode != 0) {
        // firstNode not being 0 in the continuous = false case leads to a segmentation fault
        WARN("firstNode set to 0 since continuous is false");
        this->firstNode = 0;
    }
}

const std::map<std::string, node> &EdgeListReader::getNodeMap() const {
    if (this->continuous)
        throw std::runtime_error("Input files are assumed to have continuous node ids, therefore "
                                 "no node mapping has been created.");
    return this->mapNodeIds;
}

Graph EdgeListReader::read(const std::string &path) {
    this->mapNodeIds.clear();
    MemoryMappedFile mmfile(path);
    auto it = mmfile.cbegin();
    auto end = mmfile.cend();

    bool weighted = false;
    bool checkedWeighted = false;
    Graph graph(0, weighted, directed);

    DEBUG("separator: ", this->separator);
    DEBUG("first node: ", this->firstNode);

    auto skipWhitespaceAndSeparator = [this, &it, &end]() -> void {
        while (it != end && (*it == ' ' || *it == separator))
            ++it;
    };

    // This function parses in whole words
    auto scanWord = [this, &it, &end]() -> std::string {
        std::string word = "";
        while (it != end && *it != ' ' && *it != separator && *it != '\n' && *it != '\r') {
            word += *it;
            ++it;
        }
        return word;
    };

    auto scanId = [&, this]() -> node {
        if (continuous) {
            char *past;
            auto value = strtol(it, &past, 10);
            if (past <= it)
                throw std::runtime_error("Scanning node failed. The file may be corrupt.");
            it = past;
            return value;
        } else {
            auto word = scanWord();
            auto it = mapNodeIds.find(word);
            if (it != mapNodeIds.end())
                return it->second;
            auto result = mapNodeIds.insert({word, graph.addNode()});
            if (!result.second)
                throw std::runtime_error("Error in mapping nodes");
            return result.first->second;
        }
    };

    auto scanWeight = [&it]() -> edgeweight {
        char *past;
        auto value = strtod(it, &past);
        if (past <= it)
            throw std::runtime_error("Error in parsing file - looking for weight failed");
        it = past;
        return value;
    };

    // Returns
    // 0 iff it does not point to a line ending
    // 1 iff it is a single-char ending ('\r' or '\n')
    // 2 iff it is two-char ending ("\r\n")
    auto detectLineEnding = [&it, &end]() -> size_t {
        if (*it == '\r') {
            if ((it + 1) != end && *(it + 1) == '\n')
                return 2;
            return 1;
        }
        if (*it == '\n')
            return 1;

        return 0;
    };

    // Map nodes and increase graph size
    auto mapNode = [&](node in) -> node {
        if (graph.upperNodeIdBound() <= in - firstNode)
            graph.addNodes(in - firstNode - graph.upperNodeIdBound() + 1);

        return in - firstNode;
    };

    // This function modifies the graph on input.
    auto handleEdge = [&graph](node source, node target, edgeweight weight) -> void {
        if (!graph.hasEdge(source, target)) {
            graph.addEdge(source, target, weight);
        }
    };

    while (it < end) {
        skipWhitespaceAndSeparator();
        if (detectLineEnding()) {
            // Comment line
        } else if (*it == commentPrefix[0]) {
            while (it != end && !detectLineEnding())
                ++it;
        } else {
            auto sourceId = scanId();
            if (it == end)
                throw std::runtime_error("Unexpected end of file");
            if (!(*it == ' ' || *it == separator))
                throw std::runtime_error(
                    "Error in parsing file - pointer is whitespace or separator");
            skipWhitespaceAndSeparator();

            auto targetId = scanId();
            skipWhitespaceAndSeparator();

            if (!checkedWeighted) {
                checkedWeighted = true;
                if (!detectLineEnding()) {
                    weighted = true;
                    DEBUG("Detected graph as weighted");
                    graph = GraphTools::toWeighted(graph);
                }
            }
            if (weighted) {
                auto edgeWeight = scanWeight();
                skipWhitespaceAndSeparator();

                handleEdge(mapNode(sourceId), mapNode(targetId), edgeWeight);
            } else {
                handleEdge(mapNode(sourceId), mapNode(targetId), defaultEdgeWeight);
            }
        }
        ++it;
    }

    graph.shrinkToFit();
    return graph;
}

} /* namespace NetworKit */
