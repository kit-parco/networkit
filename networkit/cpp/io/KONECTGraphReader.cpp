// no-networkit-format
/*
 * KONECTGraphReader.cpp
 *
 *  Created on: 11.05.2018
 *      Author: Roman Bange
 *
 * The KONECT format is described in detail in
 * http://konect.uni-koblenz.de/downloads/konect-handbook.pdf
 *
 */

#include <tlx/unused.hpp>

#include <networkit/auxiliary/Log.hpp>
#include <networkit/io/KONECTGraphReader.hpp>
#include <networkit/io/MemoryMappedFile.hpp>

namespace NetworKit{
    KONECTGraphReader::KONECTGraphReader(bool remapNodes, MultipleEdgesHandling handlingmethod):
    remapNodes(remapNodes), multipleEdgesHandlingMethod(handlingmethod){}

    Graph KONECTGraphReader::read(const std::string &path){
        std::string graphFormat = "";
        std::string graphType = "";
        count numberOfNodes = -1;
        count numberOfEdges = -1;

        bool directed = true;
        bool weighted = false;
        bool multiple = false;
        int valuesPerLine = 2;
        bool secondPropertyLine = false;
        std::unordered_map<node, node> nodeIdMap;

        MemoryMappedFile mmfile(path);
        auto it = mmfile.cbegin();
        auto end = mmfile.cend();

        // Returns
        // 0 iff it does not point to a line ending
        // 1 iff it is a single-char ending ('\r' or '\n')
        // 2 iff it is two-char ending ("\r\n")
        auto detectLineEnding = [&] () -> size_t {
            if (*it == '\r') {
                if ((it + 1) != end && *(it + 1) == '\n') return 2;
                return 1;
            }
            if (*it == '\n')
                return 1;

            return 0;
        };

        // Returns true if a line ending was found and skipped
        auto skipLineEnding = [&](bool required) -> bool {
            auto length = detectLineEnding();
            if (length) {
                it += length;
                return true;
            }

            if (required)
                throw std::runtime_error("No break symbol after first property line");

            return false;
        };

        // The following functions are helpers for parsing.
        auto skipWhitespace = [&] {
            while(it != end && (*it == ' ' || *it == '\t'))
                ++it;
            if(it == end){ // proper error message if file ends unexpected
                throw std::runtime_error("Unexpected end of file. Whitespace at end of file found");
            }
        };

        // This function parses in whole words
        auto scanWord = [&] () {
            std::string word = "";
            while(it != end && *it != ' ' && *it != '\t' && *it != '\n' && *it != '\r'){
                word += *it;
                ++it;
            }
            return word;
        };

        auto scanId = [&] () -> node {
            char *past;
            auto value = strtol(it, &past, 10);
            if(past <= it){
                throw std::runtime_error("Scanning node failed. The file may be corrupt.");
            }
            it = past;
            return value;
        };

        auto scanWeight = [&] () -> edgeweight{
            char *past;
            auto value = strtod(it, &past);
            if(past <= it){
                throw std::runtime_error("Scanning node failed. The file may be corrupt.");
            }
            it = past;
            return value;
        };

        //parse graph properties
        if(*it == '%'){
            ++it;
            skipWhitespace();

            // Parse graph format directed / undirected
            graphFormat = scanWord();

            if (graphFormat == "sym" || graphFormat == "bip"){
                directed = false;
            }else if(graphFormat != "asym"){
                throw std::runtime_error("Graph format not tagged properly. Format is: " + graphFormat);
            }
            skipWhitespace();
            // Parse graph weighting
            graphType = scanWord();
            //NOTE: This parser disregards the edge weight classification of "interval scale" and "ratio scale"
            if (graphType == "weighted" || graphType == "posweighted" || graphType == "signed"){
                weighted = true;
                valuesPerLine = 3;
            } else if (graphType == "positive"){ //multiple edges
                if(multipleEdgesHandlingMethod == SUM_WEIGHTS_UP){
                    weighted = true;
                }
                multiple = true; //weights will be added
            } else if (graphType == "multisigned" || graphType == "multiweighted" || graphType == "multiposweighted"){
                weighted = true;
                multiple = true;
                valuesPerLine = 3;
            } else if (graphType == "dynamic"){
                throw std::runtime_error("Dynamic graphs are not supported yet");
            } else if (graphType != "unweighted"){
                throw std::runtime_error("Graph weight not tagged properly. Weight is: " + graphType);
            }
            DEBUG("First property line read in. Format: "+graphFormat+" / Type: "+graphType);
            if (graphFormat == "bip"){
                INFO("Your graph is flagged as a bipartite one. Keep in mind that"
                 "NetworKit cannot handle this kind of graph in its special cases."
                 "It is imported as a usual undirected graph and edges might be discarded.");
            }
            if(multiple){
                DEBUG("Selected handling method for multiple edges is: "+std::to_string(multipleEdgesHandlingMethod));
            }
        } else {
            throw std::runtime_error("No graph properties line found. This line is mandatory.");
        }

        skipWhitespace();
        skipLineEnding(true);
        skipWhitespace();
        //second optional property line
        if(*it == '%'){
            ++it;
            skipWhitespace();
            numberOfEdges = scanId();
            skipWhitespace();
            numberOfNodes = scanId();
            secondPropertyLine = true;
            while(it != end && !detectLineEnding())
                ++it;
            if(it >= end){ // proper error message if file ends unexpected
                throw std::runtime_error("Unexpected end of file");
            }

            skipLineEnding(true);
            DEBUG("Second property line read in. Edges: "+std::to_string(numberOfEdges)+ " / Nodes: "+std::to_string(numberOfNodes));
            tlx::unused(numberOfEdges);
        }


        Graph graph((secondPropertyLine ? numberOfNodes : 0), weighted, directed);

        //Map nodes and increase graph size if no second property is defined
        std::function<node(node)> mapNode = [&] (node in) -> node {
            if(secondPropertyLine){
                if (in > numberOfNodes){ //if file is corrupted
                    //secondPropertyLine is made useless
                    ERROR("Given number of nodes by file does not match actual graph");
                    secondPropertyLine = false;
                    if(remapNodes){ // if remapNodes is selected true, the map has to be initalized with the existing nodes
                        nodeIdMap.reserve(numberOfNodes);
                        graph.forNodes([&](node u) {
                            nodeIdMap.insert({u,u});
                        });
                    }
                    return mapNode(in);
                }else{
                    return in - 1; //minus firstNode
                }
            }else{
                if(remapNodes){
                    auto it = nodeIdMap.find(in);
                    if(it != nodeIdMap.end())
                        return it->second;
                    auto result = nodeIdMap.insert({in, graph.addNode()});
                    assert(result.second);
                    return result.first->second;
                }else{
                    for(count i = graph.upperNodeIdBound(); i < in; i++)
                        graph.addNode();
                    return in - 1;
                }
            }
        };

        //Helper function for handling edges
        auto handleEdge = [&] (node source, node target, edgeweight weight) {
            if(!graph.hasEdge(source,target)){
                graph.addEdge(source, target, weight);
            }else if (multiple){
                switch(multipleEdgesHandlingMethod){
                    case DISCARD_EDGES: break; //Do nothing
                    case SUM_WEIGHTS_UP:
                        graph.increaseWeight(source, target, weight);
                        break;
                    case KEEP_MINIMUM_WEIGHT:
                        if(graph.weight(source,target) > weight){
                            graph.setWeight(source,target, weight);
                        }
                        break;
                    default:
                        throw std::runtime_error("Invalid multipleEdgesHandlingMethod: "
                        +std::to_string(multipleEdgesHandlingMethod));
                }
            }else{
                DEBUG("["+std::to_string(source)+"->"+std::to_string(target)+
                "] Multiple edges detected but declared as: "+graphFormat+
                " and "+graphType);
            }
        };

        while(it != end){
            skipWhitespace();
            if(skipLineEnding(false)) {
                // We ignore empty lines.
                continue;
            } else if(*it == '#' || *it == '%') {
                // Skip non-linebreak characters.
                while(it != end && !detectLineEnding())
                    ++it;
            } else {
                // Normal case parsing
                auto sourceId = scanId();
                if(it >= end){ // proper error message if file ends unexpected
                    throw std::runtime_error("Unexpected end of file");
                }
                if(!(*it == ' ' || *it == '\t')){
                    throw std::runtime_error("Target ID cannot be parsed, maybe the file is corrupt");
                }
                skipWhitespace();

                auto targetId = scanId();
                if(it >= end){ // proper error message if file ends unexpected
                    throw std::runtime_error("Unexpected end of file");
                }
                skipWhitespace();
                if (valuesPerLine > 2){
                    auto edgeWeight = scanWeight();
                    handleEdge(mapNode(sourceId), mapNode(targetId), edgeWeight);
                }else{
                    handleEdge(mapNode(sourceId), mapNode(targetId), defaultEdgeWeight);
                }
            }
            //break lines
            skipWhitespace();
            skipLineEnding(true);
        }


        graph.shrinkToFit();
        return graph;
    }
} //namespace NetworKit
