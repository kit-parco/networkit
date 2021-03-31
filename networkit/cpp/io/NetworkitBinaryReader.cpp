/*
 * NetworkitBinaryReader.cpp
 *
 * Author: Charmaine Ndolo <charmaine.ndolo@b-tu.de>
 *
 */

#include <atomic>
#include <cstring>
#include <fstream>

#include <networkit/auxiliary/Log.hpp>
#include <networkit/io/NetworkitBinaryReader.hpp>
#include <networkit/io/NetworkitBinaryGraph.hpp>
#include <networkit/io/MemoryMappedFile.hpp>

namespace NetworKit {

Graph NetworkitBinaryReader::read(const std::string& path) {
    nkbg::Header header;
    nkbg::WEIGHT_FORMAT weightFormat;

    MemoryMappedFile mmfile(path);
    auto it = mmfile.cbegin();

    auto readHeader = [&] () {
        memcpy(&header.magic, it, sizeof(uint64_t));
        it += sizeof(uint64_t);
        memcpy(&header.checksum, it, sizeof(uint64_t));
        it += sizeof(uint64_t);
        memcpy(&header.features, it, sizeof(uint64_t));
        it += sizeof(uint64_t);
        memcpy(&header.nodes, it, sizeof(uint64_t));
        it += sizeof(uint64_t);
        memcpy(&header.chunks, it, sizeof(uint64_t));
        it += sizeof(uint64_t);
        memcpy(&header.offsetBaseData, it, sizeof(uint64_t));
        it += sizeof(uint64_t);
        memcpy(&header.offsetAdjLists, it, sizeof(uint64_t));
        it += sizeof(uint64_t);
        memcpy(&header.offsetAdjTranspose, it, sizeof(uint64_t));
        it += sizeof(uint64_t);
        memcpy(&header.offsetWeightLists, it, sizeof(uint64_t));
        it += sizeof(uint64_t);
        memcpy(&header.offsetWeightTranspose, it, sizeof(uint64_t));
        it += sizeof(uint64_t);
        memcpy(&header.offsetAdjIdLists, it, sizeof(uint64_t));
        it += sizeof(uint64_t);
        memcpy(&header.offsetAdjIdTranspose, it, sizeof(uint64_t));
        it += sizeof(uint64_t);
    };

    auto checkHeader = [&] () {
        if(memcmp("nkbg002", header.magic, 8)) {
            throw std::runtime_error("Reader expected another magic value");
        } else {
            directed = (header.features & nkbg::DIR_MASK);
            weightFormat = static_cast<nkbg::WEIGHT_FORMAT>((header.features & nkbg::WGHT_MASK) >> nkbg::WGHT_SHIFT);
            indexed = (header.features & nkbg::INDEX_MASK) >> nkbg::INDEX_SHIFT;
        }
    };
    readHeader();
    checkHeader();

    nodes = header.nodes;
    DEBUG("# nodes here = ", nodes);
    chunks = header.chunks;
    DEBUG("# chunks here = ", chunks);
    if(weightFormat == nkbg::WEIGHT_FORMAT::NONE) {
        weighted = false;
    } else {
        weighted = true;
    }
    Graph G(nodes, weighted, directed);
    if(indexed) G.indexEdges();
    // Read base data.
    std::vector<uint8_t> nodeFlags;
    const char *baseIt = mmfile.cbegin() + header.offsetBaseData;
    for(uint64_t i = 0; i < nodes; i++) {
        uint8_t flag;
        memcpy(&flag, baseIt, sizeof(uint8_t));
        baseIt += sizeof(uint8_t);
        if (!(flag & nkbg::DELETED_BIT)) {
            G.removeNode(i);
        }
    }

    std::vector<uint64_t> firstVert;
    firstVert.push_back(0);
    for(uint64_t ch = 1; ch < chunks; ch++) {
        uint64_t first;
        memcpy(&first, baseIt, sizeof(uint64_t));
        baseIt += sizeof(uint64_t);
        firstVert.push_back(first);
    }
    firstVert.push_back(nodes);
    // Read adjacency lists.
    const char *adjIt = mmfile.cbegin() + header.offsetAdjLists;
    const char *transpIt = mmfile.cbegin() + header.offsetAdjTranspose;
    const char *adjWghtIt = mmfile.cbegin() + header.offsetWeightLists;
    const char *transpWghtIt = mmfile.cbegin() + header.offsetWeightTranspose;
    const char *adjIdIt = mmfile.cbegin() + header.offsetAdjIdLists;
    const char *transpIdIt = mmfile.cbegin() + header.offsetAdjIdTranspose;

    uint64_t adjListSize;
    memcpy(&adjListSize, adjIt + (chunks -1) * sizeof(uint64_t), sizeof(uint64_t));
    uint64_t transposeListSize;
    memcpy(&transposeListSize, transpIt + (chunks -1) * sizeof(uint64_t), sizeof(uint64_t));

    if(!directed) {
        assert(adjListSize == transposeListSize);
    }
    G.setEdgeCount(unsafe, adjListSize);

    std::atomic<count> selfLoops{0};
    
    auto constructGraph = [&] (uint64_t c) {
        node vertex = firstVert[c];
        uint64_t off = 0;
        uint64_t transpOff = 0;
        uint64_t wghtOff = 0;
        uint64_t transWghtOff = 0;
        uint64_t indexOff = 0;
        uint64_t transIndexOff = 0;
        if(vertex) {
            memcpy(&off, adjIt + (c-1)* sizeof(uint64_t), sizeof(uint64_t));
            memcpy(&transpOff, transpIt + (c-1)* sizeof(uint64_t), sizeof(uint64_t));
            memcpy(&wghtOff, adjWghtIt + (c-1) * sizeof(uint64_t), sizeof(uint64_t));
            memcpy(&transWghtOff, transpWghtIt + (c-1) * sizeof(uint64_t), sizeof(uint64_t));
            memcpy(&indexOff, adjIdIt + (c-1) * sizeof(uint64_t), sizeof(uint64_t));
            memcpy(&transIndexOff, transpIdIt + (c-1) * sizeof(uint64_t), sizeof(uint64_t));
        }
        off += (chunks - 1) * sizeof(uint64_t);
        transpOff += (chunks - 1) * sizeof(uint64_t);
        off += sizeof(uint64_t);
        transpOff += sizeof(uint64_t);
        wghtOff += (chunks - 1) * sizeof(uint64_t);
        transWghtOff += (chunks - 1) * sizeof(uint64_t);
        indexOff += (chunks - 1) * sizeof(uint64_t);
        transIndexOff += (chunks - 1) * sizeof(uint64_t);
        uint64_t n = firstVert[c+1] - firstVert[c];
        for (uint64_t i = 0; i < n; i++) {
            uint64_t curr = vertex+i;
            uint64_t outNbrs;
            off += nkbg::varIntDecode(reinterpret_cast<const uint8_t*>(adjIt + off), outNbrs);
            uint64_t inNbrs;
            transpOff += nkbg::varIntDecode(reinterpret_cast<const uint8_t*>(transpIt + transpOff), inNbrs);
            if(!directed) {
                G.preallocateUndirected(curr, outNbrs+inNbrs);
            } else  {
                G.preallocateDirected(curr, outNbrs, inNbrs);
            }
            //Read adjacency lists.
            for (uint64_t j = 0; j < outNbrs; j++) {
                edgeid id;
                uint64_t add;
                double weight = defaultEdgeWeight;
                off += nkbg::varIntDecode(reinterpret_cast<const uint8_t*>(adjIt + off), add);
                switch(weightFormat) {
                    case nkbg::WEIGHT_FORMAT::VARINT:
                    {
                        uint64_t unsignedWeight;
                        wghtOff += nkbg::varIntDecode(reinterpret_cast<const uint8_t*>(adjWghtIt + wghtOff), unsignedWeight);
                        weight = unsignedWeight;
                    }
                        break;
                    case nkbg::WEIGHT_FORMAT::DOUBLE:
                            memcpy(&weight, adjWghtIt + wghtOff, sizeof(double));
                            wghtOff += sizeof(double);
                        break;
                    case nkbg::WEIGHT_FORMAT::SIGNED_VARINT:
                    {
                        uint64_t unsignedWeight;
                        wghtOff += nkbg::varIntDecode(reinterpret_cast<const uint8_t*>(adjWghtIt + wghtOff), unsignedWeight);
                        weight = nkbg::zigzagDecode(unsignedWeight);
                    }
                        break;
                    case nkbg::WEIGHT_FORMAT::FLOAT:
                    {
                        float floatWeight;
                        memcpy(&floatWeight, adjWghtIt + wghtOff, sizeof(float));
                        wghtOff += sizeof(float);
                        weight = floatWeight;
                    }
                        break;
                    case nkbg::WEIGHT_FORMAT::NONE:
                        break;
                }
                if(indexed){
                    indexOff += nkbg::varIntDecode(reinterpret_cast<const uint8_t*>(adjIdIt + indexOff), id);
                }
                if(!directed) {
                    if(indexed) G.addPartialEdge(unsafe, curr, add, weight, id);
                    else G.addPartialEdge(unsafe, curr, add, weight);
                } else {
                    if(indexed) G.addPartialOutEdge(unsafe, curr, add, weight, id);
                    else G.addPartialOutEdge(unsafe, curr, add, weight);
                }
                if(curr == add) {
                    selfLoops.fetch_add(1, std::memory_order_relaxed);
                }
            }
            //Read transpose lists.
            for (uint64_t j = 0; j < inNbrs; j++) {
                uint64_t add;
                edgeid id;
                double weight = defaultEdgeWeight;
                transpOff += nkbg::varIntDecode(reinterpret_cast<const uint8_t*>(transpIt + transpOff), add);
                switch(weightFormat) {
                    case nkbg::WEIGHT_FORMAT::VARINT:
                    {
                        uint64_t unsignedWeight;
                        transWghtOff += nkbg::varIntDecode(reinterpret_cast<const uint8_t*>(transpWghtIt + transWghtOff), unsignedWeight);
                        weight = unsignedWeight;
                    }
                        break;
                    case nkbg::WEIGHT_FORMAT::DOUBLE:
                        memcpy(&weight, transpWghtIt + transWghtOff, sizeof(double));
                        transWghtOff += sizeof(double);
                        break;
                    case nkbg::WEIGHT_FORMAT::SIGNED_VARINT:
                    {
                        uint64_t unsignedWeight;
                        transWghtOff += nkbg::varIntDecode(reinterpret_cast<const uint8_t*>(transpWghtIt + transWghtOff), unsignedWeight);
                        weight = nkbg::zigzagDecode(unsignedWeight);
                    }
                        break;
                    case nkbg::WEIGHT_FORMAT::FLOAT:
                    {
                        float floatWeight;
                        memcpy(&floatWeight, transpWghtIt + transWghtOff, sizeof(float));
                        transWghtOff += sizeof(float);
                        weight = floatWeight;
                    }
                        break;
                    case nkbg::WEIGHT_FORMAT::NONE:
                        break;
                }
                if(indexed){
                    transIndexOff += nkbg::varIntDecode(reinterpret_cast<const uint8_t*>(transpIdIt + transIndexOff), id);
                }
                if(!directed) {
                    if(curr != add) {
                        if(indexed) G.addPartialEdge(unsafe, curr, add, weight, id);
                        else G.addPartialEdge(unsafe, curr, add, weight);
                    }
                } else {
                    if(indexed) G.addPartialInEdge(unsafe, curr, add, weight, id);
                    else G.addPartialInEdge(unsafe, curr, add, weight);
                }
            }
        }
    };

    // create graph
    #pragma omp parallel for
    for(omp_index c = 0; c < static_cast<omp_index>(chunks); c++) {
        constructGraph(c);
    }

    G.setNumberOfSelfLoops(unsafe, selfLoops);
    return G;
}
} /* namespace */
