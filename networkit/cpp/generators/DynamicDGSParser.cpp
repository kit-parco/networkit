// no-networkit-format
/*
 * DynamicDGSParser.cpp
 *
 *  Created on: Jun 17, 2013
 *      Author: forigem
 */

#include <fstream>

#include <networkit/auxiliary/Log.hpp>
#include <networkit/auxiliary/Parallel.hpp>
#include <networkit/generators/DynamicDGSParser.hpp>

namespace NetworKit {

DynamicDGSParser::DynamicDGSParser(const std::string &path) : graphInitialized(false) {
    dgsFile.open(path.c_str(), std::ifstream::in);
}

void DynamicDGSParser::initializeGraph() {
    if (! dgsFile.is_open()) {
        throw std::runtime_error("DGS input file could not be opened.");
    }
    else {
        DEBUG("Opened DGS file");

        std::string line;
        std::string cookie = "DGS004";

        // handle first line
        std::getline(dgsFile, line);
        if (!line.compare(0, cookie.size(), cookie)) { // compare prefix
            DEBUG("found magic cookie: DGS004");
        } else {
            DEBUG("First line: " , line);
            throw std::runtime_error("This does not seem to be a valid DGS file. Expected magic cookie 'DGS004' in first line");
        }

        // handle second line: optional name of file, number of clock ticks, total number of events
        std::getline(dgsFile, line);

        // Throw away the st0
        std::getline(dgsFile, line);

        graphInitialized = true;
    }
}

void DynamicDGSParser::generate() {
    if (graphInitialized == false) {
        throw std::runtime_error("Can not call generate() before graph was initialized.");
    }
    std::string line;
    bool breakTimeStep = false; // true if breaking from the while loop was due to a time step event

    while (std::getline(dgsFile, line)) {
        std::vector<std::string> split = Aux::StringTools::split(line);
        std::string tag = split[0];


        if (tag.compare("st") == 0 && split.size() == 2) { // clock
            Gproxy->timeStep();
            breakTimeStep = true;
            break;

        } else if (tag.compare("an") == 0 && split.size() >= 2) { // add node
            // Get the node name from the input
            std::string nodeName = split[1];
            // Add a node to a graph, mapping it to the node name inside the nodeNames map
            nodeNames[nodeName] = Gproxy->addNode();
            if (split.size() >= 4) { // DGS with ground truth

                std::string categoriesFullString = split[2]; /// Example: category="cond-mat.stat-mech, q-fin.ST"
                std::vector<std::string> categoriesFullStringSplit = Aux::StringTools::split(categoriesFullString, '"');

                std::string categoriesCommaSeparated = categoriesFullStringSplit[1]; // Example: cond-mat.stat-mech, q-fin.ST
                std::vector<std::string> categories = Aux::StringTools::split(categoriesCommaSeparated, ',');

                std::vector<std::string> currentNodeCategories(categories.begin(), categories.end());
                nodeCategories.push_back(currentNodeCategories);
                assert(!nodeCategories.empty());

                std::string dateFullString = split[3]; // Example: date="08-1997"
                std::vector<std::string> dateFullStringSplit = Aux::StringTools::split(dateFullString, '"');
                std::string date = dateFullStringSplit[1];
                nodeDates.push_back(date);

            }

        } else if (tag.compare("ae") == 0 && split.size() >= 4) { // add edge
            std::string edge_from = split[2];
            std::string edge_to = split[3];
            std::string edge_name = split[1];
            Gproxy->addEdge(nodeNames[edge_from], nodeNames[edge_to], 1.0);

        } else if (tag.compare("ce") == 0 && split.size() == 3) { // update edge. Only the "weight" attribute is supported so far

            std::string from_to_edges = split[1];
            std::vector<std::string> edgesSplit = Aux::StringTools::split(from_to_edges, '-');
            std::string edge_from = edgesSplit[0];
            std::string edge_to = edgesSplit[1];

            std::string weight = split[2];
            std::vector<std::string> weightSplit = Aux::StringTools::split(weight, '=');
            double weightValue = atoi(weightSplit[1].c_str() );

            Gproxy->setWeight(nodeNames[edge_from], nodeNames[edge_to], weightValue);

        } else if (tag.compare("dn") == 0 && split.size() == 2) {
            std::string nodeName = split[1];
            node deleteNode = nodeNames[nodeName];
            // Delete the nodes only if there are no edges connected to it
            if (Gproxy->G->degree(deleteNode) == 0) {
                Gproxy->removeNode(deleteNode);
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

            Gproxy->removeEdge(u, v);
        }

    } // end while
    if (!breakTimeStep) {
        // while loop finished because it hit the end of the file
        throw std::logic_error("reached the end of the .dgs file");
    }



}

void DynamicDGSParser::evaluateClusterings(const std::string &path, const Partition &clustering) {
    // Stores the category breakdown for each of the clusterings
    std::vector<std::unordered_map<std::string, index>> clusterMappings;
    clusterMappings.reserve(100000000); //should this be clustering.numberOfClusters()?

    // Stores the mapping between  clustring IDs (seem somewhat random
    // and not sequential) to "normalized" IDSs (0..numberOfClusterings-1)
    std::unordered_map<index, index> clusterIDMap;

    index normalizedID = 0;
    index clusterIDCounter = 0;

    //Vector of cluster sizes
    std::vector<index> clusterSizes;

    // Iterating through nodes and clusters to which they belong
    clustering.forEntries([&](node v, index c){

        // Checking if the specific cluster has been already assigned with a normalized ID
        auto gotClusterID = clusterIDMap.find (c);

        // If not, assigning it one
        if ( gotClusterID == clusterIDMap.end() ) {
            clusterIDMap[c] = clusterIDCounter;
            clusterSizes.push_back(clustering.getMembers(c).size());
            clusterIDCounter++;
        }

        // Using the normalized ID to address entries in the vector
        normalizedID = clusterIDMap[c];

        TRACE("normalizedID " , normalizedID);
        TRACE("node categories size " , nodeCategories.size());
        TRACE("c " , c);
        TRACE("v " , v);

        // Assuming 1 to 1 node correspondence (VERIFY),
        // getting the categories to which this node (a scientific paper) belongs
        std::vector<std::string> currentNodeCategories = nodeCategories[v];
        clusterMappings[normalizedID].reserve(100000);

        for (const std::string &category : currentNodeCategories) {
            clusterMappings[normalizedID].find (category);

            std::unordered_map<std::string, index>::const_iterator got = clusterMappings[normalizedID].find (category);

            if ( got == clusterMappings[normalizedID].end() ) {
                clusterMappings[normalizedID][category] = 1;
            } else {
                clusterMappings[normalizedID][category] = clusterMappings[normalizedID][category] + 1;
            }
        }
    });

    std::ofstream fs;
    fs.open (path, std::ofstream::out | std::ofstream::app);

    INFO("Clustering file output");

    for (index i=0; i<clustering.numberOfSubsets(); i++) {
        if (clusterSizes[i] >= 20) {
            fs << "cluster-" << i << '\t';
            fs << "size " << '\t' << clusterSizes[i] << '\t';
            std::vector<count> vals;
            vals.reserve(clusterMappings.size());

            for(const auto &pair : clusterMappings[i]) {
                vals.push_back(pair.second);
            }

            Aux::Parallel::sort(vals.begin(), vals.end(), std::greater<count>());

            int maxCategories = 5;
            int maxCategoriesCounter = 0;
            for (const auto &val : vals) {
                if (maxCategoriesCounter < maxCategories) {
                    fs << val << '\t';
                    maxCategoriesCounter++;
                } else {
                    break;
                }
            }

            fs << '\n';
        }
    }

    fs.close();
}

} /* namespace NetworKit */
