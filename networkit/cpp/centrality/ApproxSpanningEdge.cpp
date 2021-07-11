/*
 * ApproxSpanningEdge.hpp
 *
 *  Created on: 29.09.2019
 *     Authors: Eugenio Angriman <angrimae@hu-berlin.de>
 */

#include <cassert>
#include <cmath>
#include <omp.h>
#include <queue>

#include <networkit/centrality/ApproxSpanningEdge.hpp>
#include <networkit/components/BiconnectedComponents.hpp>

namespace NetworKit {

ApproxSpanningEdge::ApproxSpanningEdge(const Graph &G, double eps) : G(G), eps(eps) {
    if (!G.numberOfEdges())
        throw std::runtime_error("Error: graph is empty!");

    if (!G.hasEdgeIds())
        throw std::runtime_error("Error: edges not indexed, use indexEdges() before.");

    delta = 1. / static_cast<double>(G.numberOfNodes());
    visitedNodes.resize(
        omp_get_max_threads(),
        std::vector<NodeStatus>(G.upperNodeIdBound(), NodeStatus::NOT_IN_COMPONENT));
    edgeScores.resize(omp_get_max_threads(), std::vector<count>(G.upperEdgeIdBound(), 0));
    parents.resize(omp_get_max_threads(), std::vector<node>(G.upperNodeIdBound(), none));
    parentsEdgeIds.resize(omp_get_max_threads(), std::vector<edgeid>(G.upperNodeIdBound(), none));
}

std::vector<double> ApproxSpanningEdge::scores() const {
    assureFinished();
    std::vector<double> scores(edgeScores[0].begin(), edgeScores[0].end());
#pragma omp parallel for
    for (omp_index i = 0; i < static_cast<omp_index>(scores.size()); ++i)
        scores[i] /= static_cast<double>(nSamples);
    return scores;
}

void ApproxSpanningEdge::sampleUST() {
    auto &status = visitedNodes[omp_get_thread_num()];
    auto &scores = edgeScores[omp_get_thread_num()];
    auto &parent = parents[omp_get_thread_num()];
    auto &parentId = parentsEdgeIds[omp_get_thread_num()];
    auto &generator = Aux::Random::getURNG();

    for (const auto &sequence : sequences) {

        if (sequence.size() == 3) {
            // Biconnected component of size 3: we pick 2 edges uniformly at random as spanning
            // tree.
            const node u = sequence[Aux::Random::integer(2)];
            for (node v : sequence) {
                if (v != u)
                    status[v] = NodeStatus::NOT_VISITED;
            }
            G.forNeighborsOf(u, [&](node, const node v, edgeweight, const edgeid eid) {
                if (status[v] == NodeStatus::NOT_VISITED)
                    ++scores[eid];
            });

            for (node v : sequence)
                status[v] = NodeStatus::NOT_IN_COMPONENT;

            continue;
        }

        // We start building the spanning tree from the first node of the
        // sequence
        status[sequence[0]] = NodeStatus::IN_SPANNING_TREE;

        // All the remaining nodes in the components need to be visited
        std::for_each(sequence.begin() + 1, sequence.end(),
                      [&](const node u) { status[u] = NodeStatus::NOT_VISITED; });

        // Iterate over the remaining nodes to create the spanning tree
        count nodesInSpanningTree = 1;
        for (auto it = sequence.begin() + 1; it != sequence.end(); ++it) {
            const node walkStart = *it;
            if (status[walkStart] == NodeStatus::IN_SPANNING_TREE)
                continue;

            node currentNode = walkStart;
            do {
                // Get a random neighbor within the component
                node randomNeighbor = none;
                edgeid randomNeighborId = none;
                do {
                    std::tie(randomNeighbor, randomNeighborId) = G.getIthNeighborWithId(
                        currentNode, std::uniform_int_distribution<count>{0, G.degree(currentNode)
                                                                                 - 1}(generator));
                } while (status[randomNeighbor] == NodeStatus::NOT_IN_COMPONENT);

                assert(randomNeighbor != none);
                parent[currentNode] = randomNeighbor;
                parentId[currentNode] = randomNeighborId;
                currentNode = randomNeighbor;
            } while (status[currentNode] != NodeStatus::IN_SPANNING_TREE);

            const node walkEnd = currentNode;
            for (currentNode = walkStart; currentNode != walkEnd;
                 currentNode = parent[currentNode]) {
                status[currentNode] = NodeStatus::IN_SPANNING_TREE;
                ++scores[parentId[currentNode]];
                ++nodesInSpanningTree;
            }

            if (nodesInSpanningTree == sequence.size())
                break;
        }

#ifndef NDEBUG
        for (node u : sequence)
            assert(status[u] == NodeStatus::IN_SPANNING_TREE);
#endif
        for (node u : sequence)
            status[u] = NodeStatus::NOT_IN_COMPONENT;
    }
}

void ApproxSpanningEdge::run() {
    const auto m = static_cast<double>(G.numberOfEdges());
    nSamples = static_cast<count>(std::ceil(std::log(2. * m / delta) / (2. * eps * eps)));
    std::vector<node> sequence;
    sequence.reserve(G.numberOfNodes());

    BiconnectedComponents bcc(G);
    bcc.run();

    auto &toVisit = visitedNodes[0];
    std::queue<node> queue;

    for (auto curComponent : bcc.getComponents()) {
        if (curComponent.size() == 2) {
            // An edge connecting two articulation nodes is always in a ST
            edgeScores[0][G.edgeId(curComponent[0], curComponent[1])] = nSamples;
            continue;
        }

        if (curComponent.size() == 3) {
            // Biconnected components with 3 nodes can be handled easily
            sequences.push_back(curComponent);
            continue;
        }

        if (curComponent.size() > 2) {
            node source = curComponent[0];
            // We set all nodes to be visited to IN_SPANNING_TREE, and we re-set them to
            // NOT_IN_COMPONENT during the BFS. After the visit, the vector will be ready
            // for the next phase of the algorithm (we do not have to re-set the visited
            // nodes to NOT_IN_COMPONENT again).
            for (const node curComponentNode : curComponent) {
                if (G.degree(curComponentNode) > G.degree(source)) {
                    source = curComponentNode;
                }
                toVisit[curComponentNode] = NodeStatus::IN_SPANNING_TREE;
            }

            queue.push(source);
            toVisit[source] = NodeStatus::NOT_IN_COMPONENT;
            count nVisited = 0;

            do {
                const node u = queue.front();
                queue.pop();
                sequence.push_back(u);
                ++nVisited;
                G.forNeighborsOf(u, [&](const node v) {
                    if (toVisit[v] != NodeStatus::NOT_IN_COMPONENT) {
                        toVisit[v] = NodeStatus::NOT_IN_COMPONENT;
                        queue.push(v);
                    }
                });
            } while (!queue.empty());

            sequences.push_back(std::move(sequence));
            sequence.clear();
            assert(nVisited == curComponent.size());
        }
    }

#pragma omp parallel for schedule(guided)
    for (omp_index i = 0; i < static_cast<omp_index>(nSamples); ++i)
        sampleUST();

    for (count t = 1; t < edgeScores.size(); ++t) {
        const auto &scores = edgeScores[t];
#pragma omp parallel for
        for (omp_index j = 0; j < static_cast<omp_index>(scores.size()); ++j)
            edgeScores[0][j] += scores[j];
    }

    hasRun = true;
}
} // namespace NetworKit
