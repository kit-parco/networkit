// no-networkit-format
/*
 * SNAPGraphReader.cpp
 *
 *  Created on: 04.05.2018
 *      Author: Alexander van der Grinten
 */

#include <networkit/auxiliary/Log.hpp>
#include <networkit/auxiliary/StringTools.hpp>
#include <networkit/io/MemoryMappedFile.hpp>
#include <networkit/io/SNAPGraphReader.hpp>

namespace NetworKit {

SNAPGraphReader::SNAPGraphReader(bool directed, const bool& remapNodes, const count& nodeCount) :
directed(directed), nodeCount(nodeCount), remapNodes(remapNodes){}

Graph SNAPGraphReader::read(const std::string &path) {
    Graph graph(0,false,directed);

    //In the actual state this parameter has very little influence on the reader performance.
    //There can be a significant boost if it is possible to reserve space in the graph initialization
    if(nodeCount != 0 && remapNodes)
        nodeIdMap.reserve(nodeCount);

    // Maps SNAP node IDs to consecutive NetworKit node IDs.
    auto mapNode = [&] (node in) -> node {
        if (remapNodes){
            auto it = nodeIdMap.find(in);
            if(it != nodeIdMap.end())
                return it->second;
            auto result = nodeIdMap.insert({in, graph.addNode()});
            if(!result.second)
                throw std::runtime_error("Error in mapping nodes");
            return result.first->second;
        }
        for(count i = graph.upperNodeIdBound(); i < in + 1; i++)
            graph.addNode();
        return in;
    };

    // This function modifies the graph on input.
    auto handleEdge = [&] (node source, node target) {
        if(!graph.hasEdge(source, target)){
            graph.addEdge(source, target);
        }else{
            DEBUG("["+std::to_string(source)+"->"+std::to_string(target)+
                "] Multiple edges detected");
        }
    };

    MemoryMappedFile mmfile(path);
    auto it = mmfile.cbegin();
    auto end = mmfile.cend();

    // The following functions are helpers for parsing.
    auto skipWhitespace = [&] {
        while(it != end && (*it == ' ' || *it == '\t'))
            ++it;
    };

    auto scanId = [&] () -> node {
        char *past;
        auto value = strtol(it, &past, 10);
        if(past <= it)
            throw std::runtime_error("Error in parsing file - looking for nodeId failed");
        it = past;
        return value;
    };

    // This loop does the actual parsing.
    while(it != end) {
        if(it >= end)
            throw std::runtime_error("Unexpected end of file");

        skipWhitespace();

        if(it == end)
            throw std::runtime_error("Unexpected end of file");

        if(*it == '\n') {
            // We ignore empty lines.
        }else if(*it == '#') {
            // Skip non-linebreak characters.
            while(it != end && *it != '\n')
                ++it;
        }else{
            auto sourceId = scanId();
            if(it == end)
                throw std::runtime_error("Unexpected end of file");
            if(!(*it == ' ' || *it == '\t'))
                throw std::runtime_error("Error in parsing file - pointer is whitespace");

            skipWhitespace();

            auto targetId = scanId();
            skipWhitespace();

            handleEdge(mapNode(sourceId), mapNode(targetId));
        }

        if(it == end)
            throw std::runtime_error("Unexpected end of file");
        if(!(*it == '\n' || *it == '\r'))
            throw std::runtime_error("Line does not end with line break");
        ++it;
    }

    graph.shrinkToFit();
    return graph;
}

} // namespace NetworKit
