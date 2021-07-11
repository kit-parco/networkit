// no-networkit-format
/*
 * DGSReader.cpp
 *
 *  Created on: 01.06.2013
 *      Author: cls
 */

#include <iostream>
#include <map>
#include <unordered_map>

#include <networkit/auxiliary/Log.hpp>
#include <networkit/auxiliary/StringTools.hpp>
#include <networkit/dynamics/GraphEventProxy.hpp>
#include <networkit/io/DGSReader.hpp>

namespace NetworKit {

void DGSReader::read(std::string path, GraphEventProxy& Gproxy) {
    std::ifstream dgsFile;
    dgsFile.open(path.c_str(), std::ifstream::in);

    if (! dgsFile.is_open()) {
        throw std::runtime_error("DGS input file could not be opened.");
    }
    else {
        DEBUG("Opened DGS file " , path);

        std::string line;

        // handle first line
        std::getline(dgsFile, line);
        if (line.compare("DGS004") == 0) {
            DEBUG("found magic cookie: DGS004");
        } else {
            DEBUG("First line: " , line);
            throw std::runtime_error("This does not seem to be a valid DGS file. Expected magic cookie 'DGS004' in first line");
        }

        // handle second line: optional name of file, number of clock ticks, total number of events
        std::getline(dgsFile, line);
        std::vector<std::string> split = Aux::StringTools::split(line);

        std::unordered_map<std::string, node> nodeNames;
        std::vector<std::vector<std::string>> nodeCategories;
        std::vector<std::string> nodeDates;


        while (std::getline(dgsFile, line)) {
            std::vector<std::string> split = Aux::StringTools::split(line);
            std::string tag = split[0];

            if (tag.compare("st") == 0 && split.size() == 2) { // clock
                Gproxy.timeStep();

            } else if (tag.compare("an") == 0 && split.size() >= 2) { // add node
                // Example string: an cond_mat_9708143 category="cond-mat.stat-mech, q-fin.ST" date="08-1997"
                // Get the node name from the input
                std::string nodeName = split[1];
                // Add a node to a graph, mapping it to the node name inside the nodeNames map
                nodeNames[nodeName] = Gproxy.addNode();

                if (split.size() == 4) { // DGS with ground truth

                    std::string categoriesFullString = split[2]; /// Example: category="cond-mat.stat-mech, q-fin.ST"
                    std::vector<std::string> categoriesFullStringSplit = Aux::StringTools::split(categoriesFullString, '"');

                    std::string categoriesCommaSeparated = categoriesFullStringSplit[1]; // Example: cond-mat.stat-mech, q-fin.ST
                    std::vector<std::string> categories = Aux::StringTools::split(categoriesCommaSeparated, ',');

                    std::vector<std::string> currentNodeCategories(categories.begin(), categories.end());
                    nodeCategories.push_back(currentNodeCategories);

                    std::string dateFullString = split[3]; // Example: date="08-1997"
                    std::vector<std::string> dateFullStringSplit = Aux::StringTools::split(dateFullString, '"');
                    std::string date = dateFullStringSplit[1];
                    nodeDates.push_back(date);
                }
            } else if (tag.compare("ae") == 0 && split.size() >= 4) { // add edge
                std::string edge_from = split[2];
                std::string edge_to = split[3];
                std::string edge_name = split[1];
                Gproxy.addEdge(nodeNames[edge_from], nodeNames[edge_to], 1.0);

            } else if (tag.compare("ce") == 0 && split.size() == 3) { // update edge. Only the "weight" attribute is supported so far

                std::string from_to_edges = split[1];
                std::vector<std::string> edgesSplit = Aux::StringTools::split(from_to_edges, '-');
                std::string edge_from = edgesSplit[0];
                std::string edge_to = edgesSplit[1];

                std::string weight = split[2];
                std::vector<std::string> weightSplit = Aux::StringTools::split(weight, '=');
                double weightValue = atoi(weightSplit[1].c_str() );

                Gproxy.setWeight(nodeNames[edge_from], nodeNames[edge_to], weightValue);

            } else if (tag.compare("dn") == 0 && split.size() == 2) {
                std::string nodeName = split[1];
                node deleteNode = nodeNames[nodeName];
                // Delete the nodes only if there are no edges connected to it
                if (Gproxy.G->degree(deleteNode) == 0) {
                    Gproxy.removeNode(deleteNode);
                } else {
                    throw std::runtime_error("The node was not deleted, since there are edges attached to it.");
                }

            } else if (tag.compare("de") == 0 && split.size() == 2) {
                std::string from_to_edges = split[1];
                std::vector<std::string> edgesSplit = Aux::StringTools::split(from_to_edges, '-');
                std::string edge_from = edgesSplit[0];
                std::string edge_to = edgesSplit[1];
                node u = nodeNames[edge_from];
                node v = nodeNames[edge_to];
                Gproxy.removeEdge(u, v);
            }
        }

        std::cout << "nodeNames length " << nodeNames.size();
        std::map<std::string, node> ordered(nodeNames.begin(), nodeNames.end());
        for(auto it = ordered.begin(); it != ordered.end(); ++it)
            std::cout << " contents " << it->second << '\n';
    }
}

} /* namespace NetworKit */
