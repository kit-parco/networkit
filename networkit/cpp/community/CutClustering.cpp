// no-networkit-format
/*
 * Author: Michael Hamann <michael.hamann@kit.edu>
 */

#include <limits>
#include <sstream>
#include <stdexcept>

#include <networkit/auxiliary/Log.hpp>
#include <networkit/community/CutClustering.hpp>
#include <networkit/components/ConnectedComponents.hpp>
#include <networkit/flow/EdmondsKarp.hpp>
#include <networkit/graph/GraphTools.hpp>

NetworKit::CutClustering::CutClustering(const Graph& G, NetworKit::edgeweight alpha) : CommunityDetectionAlgorithm(G), alpha(alpha) { }

void NetworKit::CutClustering::run() {
    Partition result(G->upperNodeIdBound());
    result.setUpperBound(G->upperNodeIdBound());

    // Create a weighted copy of G
    auto graph = GraphTools::toWeighted(*G);

    // Augment graph by an additional node t that is connected to all other nodes
    // via an edge of weight alpha
    node t = graph.addNode();

    graph.forNodes([&](node u) {
        if (u != t) {
            graph.addEdge(u, t, alpha);
        }
    });

    // Index edges (needed by Edmonds-Karp implementation)
    graph.indexEdges();

    // sort nodes by degree, this (heuristically) reduces the number of needed cut calculations
    // bucket sort
    count n = G->numberOfNodes();
    std::vector<node> sortedNodes(n);
    {
        std::vector<index> nodePos(n + 1, 0);

        G->forNodes([&](node u) {
            ++nodePos[n - G->degree(u)];
        });

        // exclusive prefix sum
        index tmp = nodePos[0];
        index sum = tmp;
        nodePos[0] = 0;

        for (index i = 1; i < nodePos.size(); ++i) {
            tmp = nodePos[i];
            nodePos[i] = sum;
            sum += tmp;
        }

        G->forNodes([&](node u) {
            sortedNodes[nodePos[n - G->degree(u)]++] = u;
        });
    }

    for (node u : sortedNodes) {
        // the source sides have the property that they are nested, i.e. a node that
        // is already in a cluster will always produce a source side that is completely
        // contained in its cluster
        if (!result.contains(u)) {
            EdmondsKarp flowAlgo(graph, u, t);
            flowAlgo.run();
            std::vector<node> sourceSet(flowAlgo.getSourceSet());

            // all nodes in the source side form a new cluster, this cluster might absorb other clusters
            for (node v : sourceSet) {
                result[v] = u;
            }
        }
    }
    this->result = std::move(result);
    hasRun = true;
}

std::map< NetworKit::edgeweight, NetworKit::Partition > NetworKit::CutClustering::getClusterHierarchy(const NetworKit::Graph &G) {
    std::map<edgeweight, Partition> result;
    edgeweight lower = 0, upper = 2;

    // for weighted networks the upper bound of alpha is the maximum edge weight +1
    if (G.isWeighted()) {
        // FIXME this could use an OpenMP parallel max reduction
        G.forEdges([&](node, node, edgeweight weight) {
            upper = std::max(weight, upper);
        });
        upper += 1;
    }

    // for unconnected networks the lower bound for the computation are the connected components
    ConnectedComponents connComp(G);
    connComp.run();

    // construct a partition that uses a node of each cluster as representative
    Partition lowerClusters(G.upperNodeIdBound());

    // the representative of each component
    std::vector<node> componentReps(connComp.numberOfComponents(), none);

    G.parallelForNodes([&](node u) {
        // we do not care which representative is actually used so we simply take the one that wins this obvious race condition
        componentReps[connComp.componentOfNode(u)] = u;
    });

    lowerClusters.setUpperBound(G.upperNodeIdBound()); // set a safe upper bound

    // assign representatives to all nodes
    G.parallelForNodes([&](node u) {
        lowerClusters[u] = componentReps[connComp.componentOfNode(u)];
    });

    // This is the result for parameter 0
    result.insert(std::make_pair(0, lowerClusters));

    // If there is more than one connected component, the whole graph is another valid lower bound
    if (connComp.numberOfComponents() > 1) {
        node rep = GraphTools::randomNode(G);
        Partition wholeGraph(G.upperNodeIdBound(), rep);
        wholeGraph.setUpperBound(rep + 1);

        // insert the partition that contains the whole graph as cluster into the hierarchy
        result.insert(std::make_pair(-1, std::move(wholeGraph))); // wholeGraph won't be used anymore
    }

    // Construct the upper bound for the calculation. It will be inserted into the hierarchy later when the lower bound
    // for the parameter range in which this clustering is returned has been determined.
    Partition upperClusters(G.upperNodeIdBound());
    upperClusters.allToSingletons();

    clusterHierarchyRecursion(G, lower, std::move(lowerClusters), upper, std::move(upperClusters), result); // moved values won't be used anymore

    return result;
}

void NetworKit::CutClustering::clusterHierarchyRecursion(const NetworKit::Graph &G, edgeweight lower, NetworKit::Partition lowerClusters, edgeweight upper, NetworKit::Partition upperClusters, std::map< edgeweight, Partition > &result) {
    while (true) {
        edgeweight middle = -1;

        // calculate cut values and cluster sizes for all clusters
        // FIXME: it might be more efficient to calculate them only on demand, especially the cut values. However without a possibility to efficiently enumerate the nodes of a cluster this is inefficient.
        std::map<index, count> upperSizes, lowerSizes;
        std::map<index, edgeweight> upperCut, lowerCut;

        G.forEdges([&](node u, node v, edgeweight ew) {
            if (upperClusters[u] != upperClusters[v]) {
                upperCut[upperClusters[u]] += ew;
                upperCut[upperClusters[v]] += ew;
            }

            if (lowerClusters[u] != lowerClusters[v]) {
                lowerCut[lowerClusters[u]] += ew;
                lowerCut[lowerClusters[v]] += ew;
            }
        });

        G.forNodes([&](node u) {
            ++upperSizes[upperClusters[u]];
            ++lowerSizes[lowerClusters[u]];
        });

        for (auto it : lowerSizes) { // between each lower cluster and its nested upper cluster that is not itself there might be another clustering
            if (it.second == upperSizes[upperClusters[it.first]])
                continue;

            count upperSize = upperSizes[upperClusters[it.first]];
            double lowerSize = static_cast<double>(it.second);

            // find the highest possible breakpoint, first candidate: the upper cluster in which the representative of the lower cluster is.
            edgeweight lowerWeight = lowerCut[it.first];
            edgeweight upperWeight = upperCut[upperClusters[it.first]];

            edgeweight possibleBreakpoint = (upperWeight - lowerWeight) / (lowerSize - upperSize);

            // Check for all nodes in the lower cluster if they are representative of an upper cluster and check if they give a better breakpoint
            G.forNodes([&](node u) { // FIXME this is inefficient, a better way to list all nodes of a cluster would be nice
                if (lowerClusters[u] == it.first && upperClusters[u] == u && u != it.first) {
                    edgeweight tmpBreakPoint = (upperCut[u] - lowerWeight) / (lowerSize - upperSizes[u]);

                    // check temporary breakpoint, if it is better use it as new value
                    if (tmpBreakPoint > possibleBreakpoint) {
                        possibleBreakpoint = tmpBreakPoint;
                        upperWeight = upperCut[u];
                        upperSize = upperSizes[u];
                    }
                }
            });

            // check if we have actually found a new breakpoint. Use an epsilon in order to be sure that the (almost) same value is not used twice.
            if (possibleBreakpoint + std::numeric_limits<edgeweight>::epsilon() < upper) {
                // also the actual alpha value needs to be higher as otherwise the lower clustering might be calculated again
                // because the value is not high enough.
                middle = possibleBreakpoint + std::numeric_limits<edgeweight>::epsilon();
                break;
            } else if (possibleBreakpoint > upper) {
                // this shouldn't happen and indicates that there are numerical inaccuracies or an error in the implementation
                WARN("Warning: possible breakpoint ", possibleBreakpoint, " higher than the upper value: ", upper);
            }
        }

        if (middle == -1) { // no breakpoints have been found, this means that the upper bound is a tight lower bound for the upper clusters
            if (result.count(upper) == 0)
                result.insert(std::make_pair(upper, std::move(upperClusters))); // upperClusters won't be used anymore

            break;
        }

        // calculate the clustering at the calculated breakpoint
        CutClustering middleClusterer(G, middle);
        middleClusterer.run();
        Partition middleClusters(middleClusterer.getPartition()); // FIXME using a single cut clustering instance such that G is not always copied is probably faster

        INFO("Calculated clustering for alpha value ", middle);

        // As all clusterings are nested we can compare by simply comparing the number of clusters.
        count numMiddleClusters = middleClusters.numberOfSubsets();
        bool lowerIsMiddle = (lowerSizes.size() == numMiddleClusters), middleIsUpper = (numMiddleClusters == upperSizes.size());

        if (lowerIsMiddle) { // by definition the new clustering can never be the lower clustering. However this can happen because of numerical inaccuracies.
            throw std::logic_error("Error: Lower clustering is middle clustering, probably numerical inaccuracies caused this");
        } else if (middleIsUpper) { // We found the upper clustering again, this means that middle is the lowest value for which the upper clustering is returned.
            upper = middle;
            // insert the upper clustering with the determined lower bound of the parameter range
            result.insert(std::make_pair(upper, upperClusters));
        } else { // We found a new clustering between the lower and the upper clustering. Use recursion for the lower part of the interval and the loop for the upper part.
            clusterHierarchyRecursion(G, lower, std::move(lowerClusters) /* lowerClusters is overriden afterwards */, middle, middleClusters, result);
            lower = middle;
            lowerClusters = std::move(middleClusters); // middleClusters won't be used anymore
        }
    }

    assert(result.count(upper));
}


std::string NetworKit::CutClustering::toString() const {
    std::stringstream stream;

    stream << "CutClustering(" << alpha << ")";
    return stream.str();
}

