/*
 * BiasedRandomWalk.hpp
 *
 *  Created on: 03.07.2020
 *      Author: Klaus Ahrens  <ahrens@informatik.hu-berlin.de>
 *
 *  adapted and reimplemented from node2vec
 *  part of snap [https://github.com/snap-stanford/snap]
 *  Copyright (c) 2007-2019, Jure Leskovec (under BSD license)
 *
 *  see [https://arxiv.org/pdf/1607.00653v1.pdf]
 *
 */

#ifndef BIASEDRANDOMWALK_HPP
#define BIASEDRANDOMWALK_HPP

#include <memory>
#include <unordered_map>
#include <unordered_set>
#include <vector>

#include <networkit/graph/Graph.hpp>

namespace NetworKit {
namespace Embedding {

class BiasedRandomWalk {
public:
    using Walk = std::vector<node>;     // one walk
    using AllWalks = std::vector<Walk>; // n walks for all nodes: 2 dimensions in one big chunk

    /// preprocesses transition probabilities for random walks. Has to be called once before doWalks
    /// calls
    void preprocessTransitionProbs(double paramP, double paramQ);

    /// Simulates walks from every node and writes it into walks vector
    AllWalks doWalks(count walkLen, count numberOfWalks);

    BiasedRandomWalk(const Graph *graph);

private:
    const Graph *graph;
    using NeighborSet = std::unordered_set<node>;
    using NeighborMap = std::unordered_map<node, AliasSampler>;

    struct GraphData {
        // assuming contiguous node numbers 0..n-1, they serve as indices
        using Data = std::vector<NeighborMap>;
        Data data;

        GraphData(size_t nn) : data{nn, NeighborMap()} {}
    };

    std::unique_ptr<GraphData> graphData;
    std::vector<std::vector<node>> index2node;

    Walk oneWalk(node start, count walkLen);
    void preprocessNode(node t, double paramP, double paramQ);

}; // class BiasedRandomWalk
} // namespace Embedding
} // namespace NetworKit

#endif // BIASEDRANDOMWALK_HPP
