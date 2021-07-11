#include <gtest/gtest.h>

#include <networkit/generators/ErdosRenyiGenerator.hpp>
#include <networkit/graph/BFS.hpp>
#include <networkit/graph/DFS.hpp>
#include <networkit/graph/Dijkstra.hpp>
#include <networkit/graph/GraphTools.hpp>

namespace NetworKit {

class TraversalGTest : public testing::TestWithParam<std::pair<bool, bool>> {
protected:
    Graph generateRandomWeights(const Graph &G) const;
    bool weighted() const noexcept;
    bool directed() const noexcept;
};

INSTANTIATE_TEST_SUITE_P(InstantiationName, TraversalGTest,
                         testing::Values(std::make_pair(false, false), std::make_pair(true, false),
                                         std::make_pair(false, true), std::make_pair(true, true)));

Graph TraversalGTest::generateRandomWeights(const Graph &G) const {
    Graph Gw(G, true, G.isDirected());
    Gw.forEdges([&](node u, node v) { Gw.setWeight(u, v, Aux::Random::probability()); });
    return Gw;
}

bool TraversalGTest::weighted() const noexcept {
    return GetParam().first;
}

bool TraversalGTest::directed() const noexcept {
    return GetParam().second;
}

TEST_P(TraversalGTest, testBFSfrom) {
    constexpr count n = 200;
    constexpr double p = 0.15;
    std::vector<unsigned char> visited(n);
    std::vector<node> sequence;
    sequence.reserve(n);
    std::vector<std::pair<node, node>> edgeSequence;

    auto doBFS = [&](const Graph &G, const std::vector<node> &sources) {
        std::fill(visited.begin(), visited.end(), 0);
        sequence.clear();
        edgeSequence.clear();
        std::queue<node> q;

        for (node source : sources) {
            q.push(source);
            visited[source] = 1;
        }

        do {
            node u = q.front();
            q.pop();
            sequence.push_back(u);
            G.forNeighborsOf(u, [&](node v) {
                if (!visited[v]) {
                    q.push(v);
                    visited[v] = 1;
                    edgeSequence.push_back({u, v});
                }
            });
        } while (!q.empty());
    };

    std::vector<node> randNodes;
    for (node u = 0; u < n; ++u) {
        randNodes.push_back(u);
    }

    for (int seed : {1, 2, 3}) {
        Aux::Random::setSeed(seed, false);
        std::shuffle(randNodes.begin(), randNodes.end(), Aux::Random::getURNG());
        const auto G = ErdosRenyiGenerator(n, p, directed()).generate();
        for (count i = 1; i <= n; ++i) {
            std::vector<node> sources(randNodes.begin(), randNodes.begin() + i);
            doBFS(G, sources);
            count curNode = 0;
            Traversal::BFSfrom(G, sources.begin(), sources.end(),
                               [&](node u) { EXPECT_EQ(sequence[curNode++], u); });

            sources.clear();
            sources.push_back(randNodes[i - 1]);
            doBFS(G, sources);
            curNode = 0;
            Traversal::BFSEdgesFrom(G, randNodes[i - 1], [&](node u, node v, edgeweight, edgeid) {
                EXPECT_EQ(edgeSequence[curNode++], std::make_pair(u, v));
            });
        }
    }
}

TEST_P(TraversalGTest, testDFSfrom) {
    constexpr count n = 200;
    constexpr double p = 0.15;
    std::vector<unsigned char> visited(n);
    std::vector<node> sequence;
    sequence.reserve(n);
    std::vector<std::pair<node, node>> edgeSequence;

    auto doDFS = [&](const Graph &G, node source) {
        sequence.clear();
        edgeSequence.clear();
        std::fill(visited.begin(), visited.end(), 0);
        visited[source] = 1;
        std::stack<node> s;
        s.push(source);

        do {
            node u = s.top();
            s.pop();
            sequence.push_back(u);
            G.forNeighborsOf(u, [&](node v) {
                if (!visited[v]) {
                    s.push(v);
                    visited[v] = 1;
                    edgeSequence.push_back({u, v});
                }
            });
        } while (!s.empty());
    };

    for (int seed : {1, 2, 3}) {
        Aux::Random::setSeed(seed, false);
        const auto G = ErdosRenyiGenerator(n, p, directed()).generate();
        G.forNodes([&](node source) {
            doDFS(G, source);
            count curNode = 0;
            Traversal::DFSfrom(G, source, [&](node u) { EXPECT_EQ(sequence[curNode++], u); });
            curNode = 0;
            Traversal::DFSEdgesFrom(G, source, [&](node u, node v, edgeweight, edgeid) {
                EXPECT_EQ(edgeSequence[curNode++], std::make_pair(u, v));
            });
        });
    }
}

TEST_P(TraversalGTest, testDijkstraFrom) {
    constexpr count n = 200;
    constexpr double p = 0.15;

    std::vector<node> randNodes, nodes;
    for (node u = 0; u < n; ++u) {
        randNodes.push_back(u);
        nodes.push_back(u);
    }

    auto dijkstra = [&](const Graph &G, const std::vector<node> &sources) {
        struct CompareDistance {
            CompareDistance(const std::vector<edgeweight> &distance) : distance(distance) {}

            bool operator()(const node x, const node y) const noexcept {
                return distance[x] < distance[y];
            }

        private:
            const std::vector<edgeweight> &distance;
        };

        std::vector<edgeweight> distance(G.upperNodeIdBound(),
                                         std::numeric_limits<edgeweight>::max());
        tlx::d_ary_addressable_int_heap<node, 2, CompareDistance> heap{CompareDistance(distance)};
        for (const auto u : sources) {
            distance[u] = 0;
            heap.push(u);
        }

        do {
            const auto u = heap.extract_top();
            G.forNeighborsOf(u, [&](const node v, const edgeweight w) {
                if (distance[v] > distance[u] + w) {
                    distance[v] = distance[u] + w;
                    heap.update(v);
                }
            });
        } while (!heap.empty());

        return distance;
    };

    for (int seed : {1, 2, 3}) {
        Aux::Random::setSeed(seed, false);
        auto G = ErdosRenyiGenerator(n, p, directed()).generate();
        if (weighted()) {
            G = generateRandomWeights(G);
        }

        G.forNodes([&](const node u) {
            const auto distance = dijkstra(G, {u});
            std::sort(nodes.begin(), nodes.end(), [&distance](const node x, const node y) {
                return distance[x] < distance[y];
            });
            for (node i = 0; i < n - 1; ++i)
                assert(distance[nodes[i]] <= distance[nodes[i + 1]]);

            index i = 0;
            Traversal::DijkstraFrom(G, u, [&](const node u, const edgeweight w) {
                if (u != nodes[i]) {
                    EXPECT_DOUBLE_EQ(distance[u], distance[nodes[i]]);
                }
                ++i;
                EXPECT_DOUBLE_EQ(w, distance[u]);
            });
        });

        std::shuffle(randNodes.begin(), randNodes.end(), Aux::Random::getURNG());
        for (count i = 1; i <= n; ++i) {
            const auto distance =
                dijkstra(G, std::vector<node>(randNodes.begin(), randNodes.begin() + i));
            Traversal::DijkstraFrom(G, randNodes.begin(), randNodes.begin() + i,
                                    [&distance](const node u, const edgeweight d) {
                                        EXPECT_DOUBLE_EQ(d, distance[u]);
                                    });
        }
    }
}

} // namespace NetworKit
