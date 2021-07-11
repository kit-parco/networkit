// no-networkit-format
/*
 * MocnikGenerator.hpp
 *
 * Created on: July 7, 2018
 * Author: Franz-Benjamin Mocnik <mail@mocnik-science.net>
 */

#ifndef NETWORKIT_GENERATORS_MOCNIK_GENERATOR_HPP_
#define NETWORKIT_GENERATORS_MOCNIK_GENERATOR_HPP_

#include <networkit/generators/StaticGraphGenerator.hpp>

namespace NetworKit {

/**
 * @ingroup generators
 */
class MocnikGenerator final: public StaticGraphGenerator {
    // GENERAL DATA

    /**
     * Position of each node in space.  The index of the vector is also the number of
     * the node.
     */
    std::vector<std::vector<double>> nodePositions;

    // DATA FOR EACH LAYER

    /**
     * Collection of nodes.
     */
    using NodeCollection = std::vector<node>;

    /**
     * The cell array reflects how nodes are assigned to a grid. Each element of
     * the vector corresponds to one grid cell.
     */
    using CellArray = std::vector<NodeCollection>;

    /**
     * State of a layer
     */
    struct LayerState {
        CellArray a;
        int aMax;
    };

    // FUNCTIONS RELATED TO THE LAYER STATE

    /**
     * Initialize the cell array.  The second parameter determines how many grid
     * cells shall be contained in each dimension of the cell array.
     */
    void initCellArray(LayerState &s, count numberOfCellsPerDimension);

    /**
     * Get all nodes that are contained in the i-th grid cell
     */
    NodeCollection getNodes(LayerState &s, int i);

    /**
     * Add the node with number j to the corresponing grid cell
     */
    void addNode(LayerState &s, int j);

    /**
     * Determine, for a given position, the index of the corresponding grid cell
     */
    int toIndex(LayerState &s, const std::vector<double> &v);

    /**
     * Determine, for the given multi-dimensional index, the index of the
     * corresponding grid cell
     */
    int toIndex(LayerState &s, const std::vector<int> &v);

    /**
     * Determine, for a given index, the multi-dimensional index of a grid cell
     */
    std::vector<int> fromIndex(LayerState &s, int i);

    /**
     * Determine for a grid cell given by index i the grid cells of distance r
     */
    std::vector<int> boxSurface(LayerState &s, int i, int r);

    /**
     * Determine for a grid cell given by index i the grid cells within distance r
     */
    std::vector<int> boxVolume(LayerState &s, int j, double r);

    // EDGE GENERATION

    /**
     * Add the edges for n nodes and with parameter k to the graph.  Thereby, use the
     * provided relative weight.  If it is indicated that the edge is part of the
     * base layer, it is not tested whether the edge already exists inside the graph.
     */
    void addEdgesToGraph(Graph &G, count n, double k, double relativeWeight, bool baseLayer);

    count dim;
    std::vector<count> ns;
    std::vector<double> ks;
    bool weighted;
    std::vector<double> relativeWeights;

public:
    /**
     * Creates random spatial graphs according to the Mocnik model.
     *
     * Please cite the following publications, in which you will find a
     * description of the model:
     *
     * Franz-Benjamin Mocnik: "The Polynomial Volume Law of Complex Networks in
     * the Context of Local and Global Optimization", Scientific Reports 8(11274)
     * 2018. doi: 10.1038/s41598-018-29131-0
     *
     * Franz-Benjamin Mocnik, Andrew Frank: "Modelling Spatial Structures",
     * Proceedings of the 12th Conference on Spatial Information Theory (COSIT),
     * 2015, pages 44-64. doi: 10.1007/978-3-319-23374-1_3
     *
     * Improved algorithm.
     *
     * @param dim  Dimension of the space.
     * @param n  Number of nodes in the graph; or a list containing the numbers
     *     of nodes in each layer in case of a hierarchical model.
     * @param k  Density parameter, determining the ratio of edges to nodes; in
     *     case of a hierarchical model, also a list of density parameters can be
     *     provided.
     * @param weighted  Determines whether weights should be added to the edges;
     *     in case of a hierarchical model, also a list of relative weights can be
     *     provided.
     */
    MocnikGenerator(count dim, count n, double k, bool weighted=false);
    MocnikGenerator(count dim, std::vector<count> ns, double k, bool weighted = false);
    MocnikGenerator(count dim, std::vector<count> ns, std::vector<double> ks, bool weighted=false);
    MocnikGenerator(count dim, count n, double k, std::vector<double> weighted);
    MocnikGenerator(count dim, std::vector<count> ns, double k,
                    std::vector<double> weighted);
    MocnikGenerator(count dim, std::vector<count> ns, std::vector<double> ks, std::vector<double> weighted);

    Graph generate() override;
};

} /* namespace NetworKit */
#endif // NETWORKIT_GENERATORS_MOCNIK_GENERATOR_HPP_
