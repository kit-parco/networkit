// no-networkit-format
/*
 * GraphToolBinaryReader.cpp
 *
 *  Created on: 02.12.14
 *      Author: Maximilian Vogel
 */

#include <networkit/auxiliary/Enforce.hpp>
#include <networkit/auxiliary/Log.hpp>
#include <networkit/io/GraphToolBinaryReader.hpp>

namespace NetworKit {

Graph GraphToolBinaryReader::read(const std::string& path) {
    std::ifstream file(path, std::ios::binary | std::ios::in);
    Aux::enforceOpened(file);
    // if the header is ok, continue reading the file
    if(checkHeader(file)) {
        // TODO: write function to ignore the comment by skipping it.
        readComment(file);
        bool directed = getDirected(file);
        count n = getNumNodes(file);
        Graph G(n,false,directed);
        addOutNeighbours(file,n,G);
        file.close();
        G.shrinkToFit();
        return G;
    } else {
        throw std::runtime_error("File header is broken");
    }

}

bool GraphToolBinaryReader::checkHeader(std::ifstream& file) {
    uint8_t header[8];
    uint8_t reference[6] = {0xe2, 0x9b, 0xbe, 0x20, 0x67, 0x74};
    file.read((char*)header,8);
    // compare the first 6 bytes of the header
    if (!std::equal(std::begin(reference),std::end(reference),std::begin(header))) {
        return false;
    }

    if (header[6] != 0x01) {
        // version byte not valid
        return false;
    }

    // get information about endianness
    if (header[7] == 0x00) {
        this->littleEndianness = true;
    } else if (header[7] == 0x01) {
        this->littleEndianness = false;
    } else {
        // endianness byte is wrong
        return false;
    }
    return true;
}

std::string GraphToolBinaryReader::readComment(std::ifstream& file) {
    // get length of comment
    auto length = readType<uint64_t>(file,8);
    char* comment = new char[length];
    // read comment to char array
    file.read((char*)comment, length);
    std::string s(comment, length);
    delete[] comment;
    return s;
}

bool GraphToolBinaryReader::getDirected(std::ifstream& file) {
    char c;
    bool directed = false;
    file.read(&c, sizeof(char));
    if (c == 0x00) {
        directed = false;
    } else if (c == 0x01) {
        directed = true;
    } // optional else case: error handling
    return directed;
}

uint64_t GraphToolBinaryReader::getNumNodes(std::ifstream& file) {
    return readType<uint64_t>(file,8);
}

uint8_t GraphToolBinaryReader::getAdjacencyWidth(uint64_t n) {
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

std::vector<std::vector<uint64_t>> GraphToolBinaryReader::getOutNeighbours(std::ifstream& file, uint64_t numNodes) {
    std::vector<std::vector<uint64_t>> adjacencies(numNodes);
    // value of numNodes determines the size of the unsigned integer type storing the node ids
    int width = (int)getAdjacencyWidth(numNodes);
    // iterate over each node
    for(uint64_t u = 0; u < numNodes; ++u) {
        // get number of adjacencies for the current node
        uint64_t numOutNeighbours = readType<uint64_t>(file, 8);
        // iterate over adjacencies of a noede
        for(uint64_t v = 0; v < numOutNeighbours; ++v) {
            // read current adjacency and add it to the adjacency list
            uint64_t node = readType<uint64_t>(file, width);
            adjacencies[u].push_back(node);
        }
    }
    return adjacencies;
}

void GraphToolBinaryReader::addOutNeighbours(std::ifstream& file, uint64_t numNodes, Graph& G) {
    // value of numNodes determines the size of the unsigned integer type storing the node ids
    int width = (int)getAdjacencyWidth(numNodes);
    // iterate over each node
    for(uint64_t u = 0; u < numNodes; ++u) {
        // get number of adjacencies for the current node
        uint64_t numOutNeighbours = readType<uint64_t>(file, 8);
        for(uint64_t v = 0; v < numOutNeighbours; ++v) {
            // read current adjacency and add edge
            uint64_t node = readType<uint64_t>(file, width);
            G.addEdge(u,node);
        }
    }
}

} /* namespace NetworKit */
