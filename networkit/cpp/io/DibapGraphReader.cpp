/*
 * DibapGraphReader.cpp
 *
 *  Created on: Jun 12, 2013
 *      Author: Henning
 */

#ifndef NETWORKIT_WINDOWS

#include <networkit/auxiliary/Log.hpp>
#include <networkit/io/DibapGraphReader.hpp>

#include <cstdio>
#include <netinet/in.h>

// codes in file headers to distinguish type
namespace NetworKit {

enum class DibapIOType : short {
    XX = ('X' << 8) | 'X',
    GI = ('G' << 8) | 'I',
    GF = ('G' << 8) | 'F',
    HI = ('H' << 8) | 'I',
    HF = ('H' << 8) | 'F',
    P2 = ('P' << 8) | '2',
    P4 = ('P' << 8) | '4',
    AA = ('A' << 8) | 'A',
    T2 = ('T' << 8) | '2',
    TE = ('T' << 8) | 'E'
};

Graph DibapGraphReader::read(const std::string &path) {
    int n, i;
    short type;
    FILE *file = NULL; // TODO: Use std::ifstream; this function won't reliably close the handle

    int V = 0;
    int dvw = 0;
    std::vector<int> vw;
    std::vector<int> of;
    std::vector<int> to;
    int dew = 0;
    std::vector<int> ew;
    int dxy = 0;
    int numE2 = 0;

    // init, try to open file
    DEBUG("reading graph in DibaP format... ");
    file = fopen(path.c_str(), "r");
    if (file == NULL) {
        throw std::runtime_error("cannot open file ");
        return 0;
    }

    // detect file type
    n = fread(&type, sizeof(short), 1, file);
    if (n != 1 || ntohs(type) != static_cast<short>(DibapIOType::GI)) {
        throw std::runtime_error("bad file structure ");
        return 0;
    }

    // read number of vertices
    n = fread(&V, sizeof(int), 1, file);
    V = ntohl(V);
    TRACE("(V=%d %d) ", V, ", ", n);
    if (n != 1) {
        throw std::runtime_error("bad file structure ");
        return 0;
    }

    // read vertex weight dimension
    n = fread(&dvw, sizeof(int), 1, file);
    dvw = ntohl(dvw);
    TRACE("(dvw=%d %d) ", dvw, ", ", n);

    if (n != 1) {
        throw std::runtime_error("bad file structure ");
        return 0;
    }

    if (dvw > 0) {
        // read vertex weights
        vw.resize(V);
        n = fread(&vw[0], sizeof(int), V * dvw, file);
        TRACE("(vw %d) ", n);
        if (n != V) {
            throw std::runtime_error("bad file structure ");
        }
        for (i = 0; i < n; i++)
            vw[i] = ntohl(vw[i]);
    }

    of.resize(V + 1);
    n = fread(&of[0], sizeof(int), V + 1, file);
    TRACE("(of %d) ", n);
    if (n != V + 1) {
        throw std::runtime_error("bad file structure ");
    }
    for (i = 0; i < n; i++)
        of[i] = ntohl(of[i]);

    numE2 = of[V];
    to.resize(numE2);
    n = fread(&to[0], sizeof(int), numE2, file);
    TRACE("(to %d) ", n);
    if (n != numE2) {
        throw std::runtime_error("bad file structure ");
    }
    for (i = 0; i < n; i++)
        to[i] = ntohl(to[i]);

    n = fread(&dew, sizeof(int), 1, file);
    dew = ntohl(dew);
    TRACE("(dew=%d %d) ", dew, ", ", n);
    if (n != 1) {
        throw std::runtime_error("bad file structure ");
    }

    if (dew > 0) {
        int numWeights = numE2 * dew;
        ew.resize(numWeights);
        n = fread(&ew[0], sizeof(int), numWeights, file);
        TRACE("(ew %d) ", n);
        if (n != numWeights) {
            throw std::runtime_error("bad file structure ");
        }
        for (i = 0; i < n; i++)
            ew[i] = ntohl(ew[i]);
    }

    n = fread(&dxy, sizeof(int), 1, file);
    dxy = ntohl(dxy);
    TRACE("(dxy=%d %d) ", dxy, ", ", n);
    if (n != 1) {
        throw std::runtime_error("bad file structure ");
    }

    if (dxy > 0) {
        int numCoords = V * dxy;
        std::vector<float> xy;
        xy.resize(numCoords);
        n = fread(xy.data(), sizeof(float), numCoords, file);
        TRACE("(xy %d) ", n);
        if (n != numCoords) {
            throw std::runtime_error("bad file structure ");
        }

        coordinates.reserve(static_cast<size_t>(V));
        for (node u = 0; u < static_cast<node>(V); ++u) {
            coordinates.emplace_back(dxy);
            coordinates.back().copyFrom(xy.data() + (u * dxy));
        }
    }

    fclose(file);

    // fill graph: FIXME: so far without node weights, extend to weights
    Graph graph(static_cast<node>(V));

    if (dew > 0) {
        for (index v = 0; v < (count)V; ++v) {
            for (index e = of[v]; e < (index)of[v + 1]; ++e) {
                if (v <= (index)to[e]) {
                    graph.addEdge(v, to[e], ew[e]);
                }
            }
        }
    } else {
        for (index v = 0; v < (count)V; ++v) {
            for (index e = of[v]; e < (index)of[v + 1]; ++e) {
                if (v <= (index)to[e]) {
                    graph.addEdge(v, to[e]);
                }
            }
        }
    }

    graph.shrinkToFit();
    return graph;
}

} /* namespace NetworKit */

#endif // NETWORKIT_WINDOWS
