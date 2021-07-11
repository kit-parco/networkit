// no-networkit-format
/*
 * GraphLayoutAlgorithm.hpp
 *
 *  Created on: Apr 19, 2016
 *      Author: Michael Wegner
 */

#ifndef NETWORKIT_VIZ_GRAPH_LAYOUT_ALGORITHM_HPP_
#define NETWORKIT_VIZ_GRAPH_LAYOUT_ALGORITHM_HPP_

#include <fstream>
#include <vector>

#include <networkit/auxiliary/Enforce.hpp>
#include <networkit/graph/Graph.hpp>
#include <networkit/viz/Point.hpp>

namespace NetworKit {

/**
 * @ingroup viz
 * Abstract base class for algorithms that compute a layout of the Graph vertices in d-dimensional space.
 * The coordinates are stored in a d-dimensional @ref Point with type @code T \endcode.
 */
template<typename T>
class GraphLayoutAlgorithm {
public:
    GraphLayoutAlgorithm(const Graph& G, count dim) : G(&G), vertexCoordinates(std::vector<Point<T>>(G.upperNodeIdBound(), Point<T>(dim))) {}
    virtual ~GraphLayoutAlgorithm() = default;

    virtual void run() = 0;

    virtual std::vector<Point<T>> getCoordinates() const {
        return vertexCoordinates;
    }

    virtual count numEdgeCrossings() const {
        if (vertexCoordinates[0].getDimensions() == 2) {
            count numCrossings = 0;
            G->forEdges([&](node u, node v, edgeweight) {
                G->forEdges([&](node p, node q, edgeweight) {
                    if ((p == u && q == v) || (p == v && q == u)) return;
                    double m1 = (vertexCoordinates[v][1] - vertexCoordinates[u][1]) / (vertexCoordinates[v][0] - vertexCoordinates[u][0]);
                    double m2 = (vertexCoordinates[q][1] - vertexCoordinates[p][1]) / (vertexCoordinates[q][0] - vertexCoordinates[p][0]);

                    double b1 = vertexCoordinates[u][1] - m1 * vertexCoordinates[u][0];
                    double b2 = vertexCoordinates[p][1] - m1 * vertexCoordinates[p][0];
                    if (m1 != m2) {
                        double xIntersect = (b2 - b1) / (m1 - m2);
                        double minXE1 = std::min(vertexCoordinates[u][0], vertexCoordinates[v][0]);
                        double minXE2 = std::min(vertexCoordinates[p][0], vertexCoordinates[q][0]);
                        double maxXE1 = std::max(vertexCoordinates[u][0], vertexCoordinates[v][0]);
                        double maxXE2 = std::max(vertexCoordinates[p][0], vertexCoordinates[q][0]);

                        if (minXE1 <= xIntersect && minXE2 <= xIntersect && xIntersect <= maxXE1 && xIntersect <= maxXE2) {
                            numCrossings++;
                        }
                    } else if (b1 == b2) {
                        numCrossings++;
                    }
                });
            });

            numCrossings /= 2;
            return numCrossings;
        }

        return 0;
    }

    virtual bool writeGraphToGML(const std::string& filePath) {
        if (vertexCoordinates.empty() || vertexCoordinates[0].getDimensions() < 2
            || vertexCoordinates[0].getDimensions() > 3)
            return false;
        count dim = vertexCoordinates[0].getDimensions();
        std::ofstream file(filePath);
        Aux::enforceOpened(file);

        file << "graph [\n";
        if (G->isDirected()) {
            file << "  directed 1\n";
        }

        G->forNodes([&](node u) {
            file << "  node [\n";
            file << "    id " << u << "\n";
            file << "    graphics\n";
            file << "    [ x " << 50*vertexCoordinates[u][0] << "\n";
            file << "      y " << 50*vertexCoordinates[u][1] << "\n";
            if (dim == 3) {
                file << "      z " << 50*vertexCoordinates[u][2] << "\n";
            }
            file << "    ]\n";
            file << "  ]\n";
        });

        G->forEdges([&](node u, node v) {
            file << "  edge [\n";
            file << "    source "<< u << "\n";
            file << "    target "<< v << "\n";
            file << "  ]\n";
        });
        file << "]\n";

        file.close();

        return true;
    }

    virtual bool writeKinemage(const std::string& filePath) {
        if (vertexCoordinates.empty() || vertexCoordinates[0].getDimensions() != 3)
            return false;
        std::string fileName = filePath.substr(filePath.find_last_of("/"));
        std::ofstream file(filePath);
        Aux::enforceOpened(file);

        file << "@whitebackground" << std::endl;
        file << "@zoom 1.0" << std::endl;
        file << "@zslab 240" << std::endl;
        file << "@center 0 0 0" << std::endl;
        file << "@master{points}" << std::endl;
        file << "@group{" << fileName << "}" << std::endl;
        file << "@balllist {a} color= blue master={points} radius= 0.05" << std::endl;

        G->forNodes([&](node u) {
            file << "{a}" << vertexCoordinates[u][0] << " " << vertexCoordinates[u][1] << " " << vertexCoordinates[u][2] << std::endl;
        });

        // edges
        file << std::endl;
        file << "@subgroup {edges} dominant" << std::endl;
        file << "@vectorlist {edges} color= white" << std::endl;
        G->forEdges([&](node u, node v) {
            //if (u <= v) { // draw graph undirected
                file << "P " << vertexCoordinates[u][0] << " " << vertexCoordinates[u][1] << " " << vertexCoordinates[u][2] << std::endl;
                file << vertexCoordinates[v][0] << " " << vertexCoordinates[v][1] << " " << vertexCoordinates[v][2] << std::endl;
            //}
        });

        file << std::endl;

        file.close();
        return true;
    }

protected:
    const Graph* G;
    std::vector<Point<T>> vertexCoordinates;
};

} /* namespace NetworKit */

#endif // NETWORKIT_VIZ_GRAPH_LAYOUT_ALGORITHM_HPP_
