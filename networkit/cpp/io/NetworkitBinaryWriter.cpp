// no-networkit-format
/*
 * NetworkitBinaryWriter.cpp
 *
 * @author Charmaine Ndolo <charmaine.ndolo@hu-berlin.de>
 */

#include <cstring>
#include <fstream>

#include <tlx/math/clz.hpp>

#include <networkit/auxiliary/Enforce.hpp>
#include <networkit/auxiliary/Log.hpp>
#include <networkit/io/NetworkitBinaryGraph.hpp>
#include <networkit/io/NetworkitBinaryWriter.hpp>

namespace NetworKit {

NetworkitBinaryWriter::NetworkitBinaryWriter(uint64_t chunks, NetworkitBinaryWeights weightsType): chunks(chunks), weightsType(weightsType) {}

void NetworkitBinaryWriter::write(const Graph &G, const std::string &path) {

    std::ofstream outfile(path, std::ios::binary);
    Aux::enforceOpened(outfile);
    nkbg::WEIGHT_FORMAT weightFormat;

    auto detectWeightsType = [&] () -> nkbg::WEIGHT_FORMAT {
        if(!G.isWeighted())
            return nkbg::WEIGHT_FORMAT::NONE;

        bool isUnsigned = true;
        bool fitsIntoInt64 = true;
        bool fitsIntoFloat = true;
        G.forEdges([&](node, node, edgeweight w) {
            if(w < 0)
                isUnsigned = false;
            if(w != static_cast<int64_t>(w))
                fitsIntoInt64 = false;
            if(w != static_cast<float>(w))
                fitsIntoFloat = false;
        });
        if(fitsIntoInt64) {
            if(isUnsigned) {
                return nkbg::WEIGHT_FORMAT::VARINT;
            } else {
                return nkbg::WEIGHT_FORMAT::SIGNED_VARINT;
            }
        } else if(fitsIntoFloat) {
            return nkbg::WEIGHT_FORMAT::FLOAT;
        } else {
            return nkbg::WEIGHT_FORMAT::DOUBLE;
        }
    };

    switch(weightsType) {
        case NetworkitBinaryWeights::none:
            weightFormat = nkbg::WEIGHT_FORMAT::NONE;
            break;
        case NetworkitBinaryWeights::autoDetect:
            weightFormat = detectWeightsType();
            break;
        case NetworkitBinaryWeights::unsignedFormat:
            weightFormat = nkbg::WEIGHT_FORMAT::VARINT;
            break;
        case NetworkitBinaryWeights::signedFormat:
            weightFormat = nkbg::WEIGHT_FORMAT::SIGNED_VARINT;
            break;
        case NetworkitBinaryWeights::doubleFormat:
            weightFormat = nkbg::WEIGHT_FORMAT::DOUBLE;
            break;
        case NetworkitBinaryWeights::floatFormat:
            weightFormat = nkbg::WEIGHT_FORMAT::FLOAT;
            break;
    }

    nkbg::Header header;
    uint64_t adjListSize = 0;
    uint64_t adjTransposeSize = 0;

    auto setFeatures = [&] () {
        header.features =
            (G.isDirected() & nkbg::DIR_MASK)
            | ((static_cast<uint64_t>(weightFormat) << nkbg::WGHT_SHIFT) & nkbg::WGHT_MASK);
    };

    auto writeHeader = [&] () {
        outfile.write(header.magic, 8);
        outfile.write(reinterpret_cast<char*>(&header.checksum), sizeof(uint64_t));
        outfile.write(reinterpret_cast<char*>(&header.features), sizeof(uint64_t));
        outfile.write(reinterpret_cast<char*>(&header.nodes), sizeof(uint64_t));
        outfile.write(reinterpret_cast<char*>(&header.chunks), sizeof(uint64_t));
        outfile.write(reinterpret_cast<char*>(&header.offsetBaseData), sizeof(uint64_t));
        outfile.write(reinterpret_cast<char*>(&header.offsetAdjLists), sizeof(uint64_t));
        outfile.write(reinterpret_cast<char*>(&header.offsetAdjTranspose), sizeof(uint64_t));
        outfile.write(reinterpret_cast<char*>(&header.offsetWeightLists), sizeof(uint64_t));
        outfile.write(reinterpret_cast<char*>(&header.offsetWeightTranspose), sizeof(uint64_t));
    };

    auto writeWeightsToFile = [&] (edgeweight w) {
        uint8_t tmp [10];
        uint64_t weightSize;
        switch(weightFormat) {
            case nkbg::WEIGHT_FORMAT::VARINT:
                weightSize = nkbg::varIntEncode(w, tmp);
                outfile.write(reinterpret_cast<char*>(tmp), weightSize);
                break;
            case nkbg::WEIGHT_FORMAT::DOUBLE:
                outfile.write(reinterpret_cast<char*>(&w), sizeof(double));
                break;
            case nkbg::WEIGHT_FORMAT::SIGNED_VARINT:
            {
                uint64_t weight = nkbg::zigzagEncode(w);
                weightSize = nkbg::varIntEncode(weight,tmp);
                outfile.write(reinterpret_cast<char*>(tmp), weightSize);
            }
                break;
            case nkbg::WEIGHT_FORMAT::FLOAT:
            {
                float weight = w;
                outfile.write(reinterpret_cast<char*>(&weight), sizeof(float));
            }
                break;
            case nkbg::WEIGHT_FORMAT::NONE:
                break; // avoid compiler warning
        }
    };

    count nodes = G.numberOfNodes();
    if (nodes < chunks) {
        chunks = nodes;
        INFO("reducing chunks to ", chunks, " chunks");
    }

    // Compute first vertex of each chunk.
    std::vector<uint64_t> firstInChunk;
    firstInChunk.push_back(0);
    uint64_t firstNode = 0;
    for(uint64_t c = 1; c < chunks; c++) {
        firstNode += (nodes/chunks);
        firstInChunk.push_back(firstNode);
    }
    firstInChunk.push_back(nodes);

    // Compute encoded size of arrays and store in vector.
    uint64_t adjSize = 0;
    uint64_t transpSize = 0;
    uint64_t adjWeightSize = 0;
    uint64_t transpWeightSize = 0;
    std::vector<uint64_t> nrOutNbrs;
    std::vector<uint64_t> nrInNbrs;
    std::vector<size_t> adjOffsets; //Prefix sum of size encoded adj arrays
    std::vector<size_t> transpOffsets; //Prefix sum of encoded transposed adj arrays
    std::vector<size_t> adjWghtOffsets; //Prefix sum of size encoded adj weights
    std::vector<size_t> transpWghtOffsets;//Prefix sum of encoded transposed weights

    auto computeWeightsOffsets = [&] (edgeweight w) {
        uint64_t size = 0;
        uint8_t tmp [10];
        switch(weightFormat) {
            case nkbg::WEIGHT_FORMAT::VARINT:
                size = nkbg::varIntEncode(w,tmp);
                break;
            case nkbg::WEIGHT_FORMAT::DOUBLE:
                size = sizeof(double);
                break;
            case nkbg::WEIGHT_FORMAT::SIGNED_VARINT:
            {
                uint64_t weight = nkbg::zigzagEncode(w);
                size += nkbg::varIntEncode(weight,tmp);
            }
                break;
            case nkbg::WEIGHT_FORMAT::FLOAT:
                size += sizeof(float);
                break;
            case nkbg::WEIGHT_FORMAT::NONE:
                break; // avoid compiler warning
        }

        return size;
    };

    for(uint64_t c = 0; c < chunks; c++) {
        for(uint64_t n = firstInChunk[c]; n < firstInChunk[c+1]; n++) {
            uint64_t outNbrs = 0;
            uint64_t inNbrs = 0;
            uint8_t tmp [10];
            if(!G.isDirected()){
                G.forNeighborsOf(n,[&](node v, edgeweight w) {
                    if(v <= n) {
                        outNbrs++;
                        adjSize += nkbg::varIntEncode(v, tmp);
                        adjWeightSize += computeWeightsOffsets(w);
                    }
                    if (v >= n) {
                        inNbrs++;
                        transpSize += nkbg::varIntEncode(v, tmp);
                        transpWeightSize += computeWeightsOffsets(w);
                    }
                });
            } else {
                G.forNeighborsOf(n, [&] (node v, edgeweight w) {
                    outNbrs++;
                    adjSize += nkbg::varIntEncode(v, tmp);
                    adjWeightSize += computeWeightsOffsets(w);
                });
                G.forInNeighborsOf(n, [&] (node v, edgeweight w) {
                    inNbrs++;
                    transpSize += nkbg::varIntEncode(v, tmp);
                    transpWeightSize += computeWeightsOffsets(w);
                });
            }
            adjListSize += outNbrs;
            adjSize += nkbg::varIntEncode(outNbrs, tmp);
            nrOutNbrs.push_back(outNbrs);

            adjTransposeSize += inNbrs;
            transpSize += nkbg::varIntEncode(inNbrs, tmp);
            nrInNbrs.push_back(inNbrs);
        }
        adjOffsets.push_back(adjSize);
        transpOffsets.push_back(transpSize);
        adjWghtOffsets.push_back(adjWeightSize);
        transpWghtOffsets.push_back(transpWeightSize);
    }
    // Write header.
    strncpy(header.magic,"nkbg002",8);
    header.checksum = 0;
    setFeatures();
    header.nodes = nodes;
    header.chunks = chunks;
    header.offsetBaseData = sizeof(nkbg::Header);
    header.offsetAdjLists = header.offsetBaseData
            + nodes * sizeof(uint8_t) // nodeFlags.
            + (chunks - 1) * sizeof(uint64_t); // firstVertex.
    header.offsetAdjTranspose = header.offsetAdjLists
            + (chunks - 1) * sizeof(uint64_t) // adjOffsets
            + sizeof(uint64_t) // adjListSize
            + adjOffsets.back(); // Size of data
    if(weightFormat != nkbg::WEIGHT_FORMAT::NONE) {
            header.offsetWeightLists = header.offsetAdjTranspose
            + (chunks - 1) * sizeof(uint64_t) //transpOffsets
            + sizeof(uint64_t) //adjTranspSize
            + transpOffsets.back(); //Size of data
        header.offsetWeightTranspose = header.offsetWeightLists
            + (chunks - 1) * sizeof(uint64_t)
            + adjWghtOffsets.back();
    } else {
        header.offsetWeightLists = 0;
        header.offsetWeightTranspose = 0;
    }
    writeHeader();
    // Write base data.
    G.forNodes([&] (node u) {
        uint8_t nodeFlag = 0;
        if (G.hasNode(u)) {
            nodeFlag |= nkbg::DELETED_BIT;
        }
        outfile.write(reinterpret_cast<char*>(&nodeFlag), sizeof(uint8_t));
    });

    assert(!firstInChunk[0]);
    for (uint64_t c = 1; c < chunks; c++) {
        outfile.write(reinterpret_cast<char*>(&firstInChunk[c]), sizeof(uint64_t));
    }

    // Write adjacency data.
    for (uint64_t c = 1; c < chunks; c++) {
        outfile.write(reinterpret_cast<char*>(&adjOffsets[c-1]), sizeof(uint64_t));
    }
    // Write size of list
    outfile.write(reinterpret_cast<char*>(&adjListSize), sizeof(uint64_t));
    G.forNodes([&](node u) {
        uint8_t tmp [10];
        uint64_t nbrsSize = nkbg::varIntEncode(nrOutNbrs[u], tmp);
        outfile.write(reinterpret_cast<char*>(tmp), nbrsSize);
        G.forNeighborsOf(u,[&](node v){
            uint64_t nodeSize;
            if(!G.isDirected()) {
                if(v <= u) {
                    nodeSize = nkbg::varIntEncode(v, tmp);
                    outfile.write(reinterpret_cast<char*>(tmp), nodeSize);
                }
            } else {
                nodeSize = nkbg::varIntEncode(v, tmp);
                outfile.write(reinterpret_cast<char*>(tmp), nodeSize);
            }
        });
    });

    // Write transpose data.
    for (uint64_t c = 1; c < chunks; c++) {
        outfile.write(reinterpret_cast<char*>(&transpOffsets[c-1]), sizeof(uint64_t));
    }
    // Write size of transpose list.
    outfile.write(reinterpret_cast<char*>(&adjTransposeSize), sizeof(uint64_t));
    G.forNodes([&](node u) {
        uint8_t tmp [10];
        uint64_t nbrsSize = nkbg::varIntEncode(nrInNbrs[u], tmp);
        outfile.write(reinterpret_cast<char*>(tmp), nbrsSize);
        uint64_t nodeSize;
        if(!G.isDirected()) {
            G.forNeighborsOf(u,[&](node v){
                if(v >= u) {
                    nodeSize = nkbg::varIntEncode(v, tmp);
                    outfile.write(reinterpret_cast<char*>(tmp), nodeSize);
                }
            });
        } else {
            G.forInNeighborsOf(u,[&](node v){
                nodeSize = nkbg::varIntEncode(v, tmp);
                outfile.write(reinterpret_cast<char*>(tmp), nodeSize);
            });
        }
    });

    // Write adj weights.
    for (uint64_t c = 1; c < chunks; c++) {
        outfile.write(reinterpret_cast<char*>(&adjWghtOffsets[c-1]), sizeof(uint64_t));
    }
    G.forNodes([&](node u) {
        G.forNeighborsOf(u,[&](node v,edgeweight w){
            if(!G.isDirected()) {
                if(v <= u) {
                    writeWeightsToFile(w);
                }
            } else {
                writeWeightsToFile(w);
            }
        });
    });

    // Write transpose weights.
    for (uint64_t c = 1; c < chunks; c++) {
        outfile.write(reinterpret_cast<char*>(&transpWghtOffsets[c-1]), sizeof(uint64_t));
    }
    G.forNodes([&](node u) {
        if(!G.isDirected()) {
            G.forNeighborsOf(u,[&](node v, edgeweight w){
                if(v >= u) {
                    writeWeightsToFile(w);
                }
            });
        } else {
            G.forInNeighborsOf(u, [&](node, edgeweight w){
                writeWeightsToFile(w);
            });
        }
    });

    INFO("Written graph to ", path);
}
} /*namespace */
