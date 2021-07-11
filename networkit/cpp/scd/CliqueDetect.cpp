#include <unordered_map>

#include <networkit/auxiliary/IncrementalUniformRandomSelector.hpp>
#include <networkit/clique/MaximalCliques.hpp>
#include <networkit/graph/GraphTools.hpp>
#include <networkit/scd/CliqueDetect.hpp>

namespace NetworKit {

CliqueDetect::CliqueDetect(const Graph &g) : SelectiveCommunityDetector(g) {
    if (g.numberOfSelfLoops() > 0)
        throw std::runtime_error("CliqueDetect works only with simple graphs.");
    if (g.isDirected())
        throw std::runtime_error("CliqueDetect work only with undirected graphs.");
}

std::set<node> CliqueDetect::expandOneCommunity(node s) {
    // 1. get the maximum clique in neighbors(s) + s
    std::set<node> result;
    result.insert(s);

    std::vector<node> sn(g->neighborRange(s).begin(), g->neighborRange(s).end());
    std::vector<edgeweight> neighborWeights;
    if (g->isWeighted()) {
        neighborWeights.reserve(sn.size());

        for (auto it : g->weightNeighborRange(s)) {
            neighborWeights.push_back(it.second);
        }
    }

    if (!sn.empty()) {
        for (node u : getMaximumWeightClique(sn, neighborWeights)) {
            result.insert(sn[u]);
        }
    }

    return result;
}

std::set<node> CliqueDetect::expandOneCommunity(const std::set<node> &seeds) {
    std::set<node> result(seeds);

    if (!seeds.empty()) {
        // Find neighbors of the seeds that are neighbors of all seed nodes
        // First, candidates are neighbors of the first seed node
        std::unordered_map<node, std::pair<count, edgeweight>> sn;
        g->forNeighborsOf(*seeds.begin(), [&](node v, edgeweight weight) {
            if (seeds.find(v) == seeds.end()) {
                sn.insert({v, {1u, weight}});
            }
        });

        // Count for each neighbor to how many other seed nodes it is adjacent
        for (auto it = ++seeds.begin(); it != seeds.end(); ++it) {
            g->forNeighborsOf(*it, [&](node v, edgeweight weight) {
                auto it = sn.find(v);
                if (it != sn.end()) {
                    ++it->second.first;
                    it->second.second += weight;
                }
            });
        }

        // Take those neighbors that are adjacent to all seed nodes
        std::vector<node> neighbors;
        std::vector<edgeweight> neighborWeights;
        neighbors.reserve(sn.size());
        if (g->isWeighted()) {
            neighborWeights.reserve(sn.size());
        }

        for (auto it : sn) {
            assert(it.second.first <= seeds.size() && "Graph has probably multi-edges!");
            if (it.second.first == seeds.size()) {
                neighbors.push_back(it.first);
                if (g->isWeighted()) {
                    neighborWeights.push_back(it.second.second);
                }
            }
        }

        if (!neighbors.empty()) {
            const std::vector<node> clique = getMaximumWeightClique(neighbors, neighborWeights);

            for (node u : clique) {
                result.insert(neighbors[u]);
            }
        }
    }

    return result;
}

std::vector<node>
CliqueDetect::getMaximumWeightClique(const std::vector<node> &nodes,
                                     const std::vector<edgeweight> &seedToNodeWeight) const {

    Graph s = GraphTools::subgraphFromNodes(*g, nodes.begin(), nodes.end(), true);
    std::vector<node> maxClique;

    Aux::IncrementalUniformRandomSelector selector;
    if (!g->isWeighted()) {
        // Select a maximum clique uniformly at random
        MaximalCliques mc(s, [&](const std::vector<node> &clique) {
            if (clique.size() < maxClique.size())
                return;
            if (clique.size() > maxClique.size()) {
                maxClique = clique;
                selector.reset();
            } else {
                if (selector.addElement()) {
                    maxClique = clique;
                }
            }
        });
        mc.run();
    } else {
        // Select a maximum weight clique uniformly at random
        double maxWeight = 0.0;
        assert(s.numberOfNodes() == s.upperNodeIdBound());
        std::vector<bool> inClique(s.upperNodeIdBound(), false);

        // clique with the greatest sum of its edgeweights
        MaximalCliques mc(s, [&](const std::vector<node> &clique) {
            double cliqueWeight = 0.0;

            for (node u : clique) {
                inClique[u] = true;
            }

            for (node u : clique) {
                for (auto neighborIt : s.weightNeighborRange(u)) {
                    cliqueWeight += inClique[neighborIt.first] * neighborIt.second;
                }

                // Add the weight from u to the seed node(s)
                cliqueWeight += seedToNodeWeight[u];

                // Reset inClique[u] to avoid counting edges twice
                inClique[u] = false;
            }

            if (cliqueWeight > maxWeight) {
                maxWeight = cliqueWeight;
                maxClique = clique;
                selector.reset();
            } else if (cliqueWeight == maxWeight) {
                if (selector.addElement()) {
                    maxClique = clique;
                }
            }
        });

        mc.run();
    }

    return maxClique;
}
} /* namespace NetworKit */
