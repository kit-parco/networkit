/*
 * ApproxElectricalCloseness.cpp
 *
 *  Created on: 17.10.2019
 *     Authors: Eugenio Angriman <angrimae@hu-berlin.de>
 *              Alexander van der Grinten <avdgrinten@hu-berlin.de>
 */

#include <numeric>
#include <omp.h>
#include <queue>

#include <networkit/algebraic/CSRMatrix.hpp>
#include <networkit/auxiliary/Random.hpp>
#include <networkit/centrality/ApproxElectricalCloseness.hpp>
#include <networkit/numerics/ConjugateGradient.hpp>
#include <networkit/numerics/LAMG/Lamg.hpp>
#include <networkit/numerics/Preconditioner/DiagonalPreconditioner.hpp>

namespace NetworKit {

ApproxElectricalCloseness::ApproxElectricalCloseness(const Graph &G, double epsilon, double kappa)
    : Centrality(G), epsilon(epsilon), delta(1.0 / static_cast<double>(G.numberOfNodes())),
      kappa(kappa), bccPtr(new BiconnectedComponents(G)) {

    if (G.isDirected())
        throw std::runtime_error("Error: the input graph must be undirected.");

    if (G.isWeighted())
        throw std::runtime_error("Error: the input graph must be unweighted.");

    if (G.numberOfNodes() < 2)
        throw std::runtime_error("Error: the graph should have at leasts two vertices");

    if (G.upperNodeIdBound() != G.numberOfNodes())
        throw std::runtime_error("Error, graph is not compact. Use getCompactedGraph in GraphTools "
                                 "before using this algorithm.");

    const count n = G.upperNodeIdBound(), threads = omp_get_max_threads();
    statusGlobal.resize(threads, std::vector<NodeStatus>(n, NodeStatus::NOT_VISITED));
    parentGlobal.resize(threads, std::vector<node>(n, none));
    approxEffResistanceGlobal.resize(threads, std::vector<int64_t>(n));
    scoreData.resize(n);
    diagonal.resize(n);
    tVisitGlobal.resize(threads, std::vector<count>(n));
    tFinishGlobal.resize(threads, std::vector<count>(n));
    generators.reserve(threads);
    for (omp_index i = 0; i < threads; ++i)
        generators.emplace_back(Aux::Random::getURNG());

    degDist.reserve(G.upperNodeIdBound());
    G.forNodes([&](node u) { degDist.emplace_back(0, G.degree(u) - 1); });

    ustChildPtrGlobal.resize(threads, std::vector<node>(n));
    ustSiblingPtrGlobal.resize(threads, std::vector<node>(n));
    bfsParent.resize(n, none);
}

count ApproxElectricalCloseness::computeNumberOfUSTs() const {
    return rootEcc * rootEcc
           * static_cast<count>(
               std::ceil(std::log(2.0 * static_cast<double>(G.numberOfEdges()) / delta)
                         / (2.0 * epsilon * epsilon * (1. - kappa) * (1. - kappa))));
}

node ApproxElectricalCloseness::approxMinEccNode() {
    auto &status = statusGlobal[0];
    std::vector<count> distance(G.upperNodeIdBound());
    std::vector<count> eccLowerBound(G.upperNodeIdBound());

    auto maxDegreeNode = [&]() -> node {
        node maxDegNode = 0;
        count maxDeg = 0;
        G.forNodes([&](const node u) {
            const auto degU = G.degree(u);
            if (degU > maxDeg) {
                maxDeg = degU;
                maxDegNode = u;
            }
        });
        return maxDegNode;
    };

    auto doBFS = [&](const node source) -> node {
        std::queue<node> q;
        q.push(source);
        status[source] = NodeStatus::VISITED;
        distance[source] = 0;
        node farthest = 0;

        do {
            const node u = q.front();
            q.pop();
            eccLowerBound[u] = std::max(eccLowerBound[u], distance[u]);
            farthest = u;

            G.forNeighborsOf(u, [&](const node v) {
                if (status[v] == NodeStatus::NOT_VISITED) {
                    q.push(v);
                    status[v] = NodeStatus::VISITED;
                    distance[v] = distance[u] + 1;
                }
            });

        } while (!q.empty());

        std::fill(status.begin(), status.end(), NodeStatus::NOT_VISITED);
        return farthest;
    };

    node source = maxDegreeNode();

    for (uint32_t i = 0; i < sweeps; ++i)
        source = doBFS(source);

    // Return node with minimum ecc lower bound
    return static_cast<node>(std::min_element(eccLowerBound.begin(), eccLowerBound.end())
                             - eccLowerBound.begin());
}

void ApproxElectricalCloseness::computeNodeSequence() {
    // We use thread 0's vector
    auto &status = statusGlobal[0];

    // Compute the biconnected components
    auto &bcc = *bccPtr;
    bcc.run();
    auto components = bcc.getComponents();

    std::queue<node> queue;
    std::vector<node> curSequence;

    for (const auto &curComponent : components) {
        // Biconnected components with 2 vertices can be handled trivially.
        if (curComponent.size() == 2) {
            sequences.push_back({curComponent[0], curComponent[1]});
            continue;
        }

        // We take the node with highest degree in the component as source
        node source = curComponent[0];
        for (node u : curComponent) {
            source = G.degree(u) > G.degree(source) ? u : source;
            status[u] = NodeStatus::VISITED;
        }

        // Start a BFS from source to determine the order of the nodes.
        // Later, we need to re-explore the graph again, so in the beginning we mark all nodes to be
        // visited as VISITED, and we set them as NOT_VISITED during the BFS. This avoids a call to
        // std::fill.
        queue.push(source);
        status[source] = NodeStatus::NOT_VISITED;

        do {
            const node u = queue.front();
            queue.pop();
            curSequence.push_back(u);
            G.forNeighborsOf(u, [&](node v) {
                if (status[v] == NodeStatus::VISITED) {
                    status[v] = NodeStatus::NOT_VISITED;
                    queue.push(v);
                }
            });
        } while (!queue.empty());

        sequences.push_back(std::move(curSequence));
        curSequence.clear();
    }

    // Set the root to the highest degree node within the largest biconnected component
    root = approxMinEccNode();

    biAnchor.resize(bcc.numberOfComponents(), none);
    biParent.resize(bcc.numberOfComponents(), none);

#ifdef NETWORKIT_SANITY_CHECKS
    G.forNodes([&](node u) { assert(status[u] == NodeStatus::NOT_VISITED); });
#endif

    // Topological order of biconnected components: tree of biconnected components starting from the
    // root's biconnected component. If the root is in multiple biconnected components, we take one
    // of them arbitrarily select one of them.
    std::queue<std::pair<node, index>> q;
    const auto &rootComps = bcc.getComponentsOfNode(root);
    q.push({root, *(rootComps.begin())});

    topOrder.reserve(bcc.numberOfComponents());
    topOrder.insert(topOrder.begin(), rootComps.begin(), rootComps.end());

    std::vector<count> distance(G.upperNodeIdBound());
    status[root] = NodeStatus::VISITED;
    rootEcc = 0;

    do {
        const auto front = q.front();
        q.pop();
        G.forNeighborsOf(front.first, [&](const node v) {
            if (status[v] == NodeStatus::NOT_VISITED) {
                distance[v] = distance[front.first] + 1;
                rootEcc = std::max(rootEcc, distance[v]);
                const auto &vComps = bcc.getComponentsOfNode(v);
                for (const node vComponentIndex : vComps) {
                    // Check if a new biconnected components has been found.
                    if (vComponentIndex != front.second && biAnchor[vComponentIndex] == none
                        && rootComps.find(vComponentIndex) == rootComps.end()) {
                        // The anchor cannot be the root, because the anchor does not have a parent.
                        // We handle biAnchor = none cases later.
                        biAnchor[vComponentIndex] = (v == root) ? none : v;
                        biParent[vComponentIndex] = front.second;
                        topOrder.push_back(vComponentIndex);
                    }
                }
                q.push({v, *(vComps.begin())});
                status[v] = NodeStatus::VISITED;
            }
        });
    } while (!q.empty());

#ifdef NETWORKIT_SANITY_CHECKS
    G.forNodes([&](node u) { assert(status[u] == NodeStatus::VISITED); });
    assert(topOrder.size() == bcc.numberOfComponents());
    assert(std::unordered_set<index>(topOrder.begin(), topOrder.end()).size()
           == bcc.numberOfComponents());
#endif
}

void ApproxElectricalCloseness::computeBFSTree() {
    // Using thread 0's vector
    auto &status = statusGlobal[0];
    std::fill(status.begin(), status.end(), NodeStatus::NOT_VISITED);

    std::queue<node> queue;
    queue.push(root);
    status[root] = NodeStatus::VISITED;

    do {
        const node currentNode = queue.front();
        queue.pop();
        G.forNeighborsOf(currentNode, [&](const node v) {
            if (status[v] == NodeStatus::NOT_VISITED) {
                status[v] = NodeStatus::VISITED;
                queue.push(v);
                bfsParent[v] = currentNode;
            }
        });
    } while (!queue.empty());

#ifdef NETWORKIT_SANITY_CHECKS
    checkBFSTree();
#endif
}

void ApproxElectricalCloseness::sampleUST() {
    // Getting thread-local vectors
    auto &status = statusGlobal[omp_get_thread_num()];
    auto &parent = parentGlobal[omp_get_thread_num()];
    auto &childPtr = ustChildPtrGlobal[omp_get_thread_num()];
    auto &siblingPtr = ustSiblingPtrGlobal[omp_get_thread_num()];
    std::fill(status.begin(), status.end(), NodeStatus::NOT_IN_COMPONENT);
    std::fill(parent.begin(), parent.end(), none);
    std::fill(childPtr.begin(), childPtr.end(), none);
    std::fill(siblingPtr.begin(), siblingPtr.end(), none);

    auto &generator = generators[omp_get_thread_num()];

    // Iterate over the biconnected components in their topological order.
    for (const auto componentIndex : topOrder) {
        // Current component, sorted by vertex degree.
        const auto &sequence = sequences[componentIndex];
        auto curAnchor = biAnchor[componentIndex];

        // Finds the parent of the current anchor node i.e., the anchor's neighbor that is in the
        // parent component.
        auto updateParentOfAnchor = [&]() -> void {
            for (const node v : G.neighborRange(curAnchor)) {
                const auto &vComps = bccPtr->getComponentsOfNode(v);
                if (vComps.find(biParent[componentIndex]) != vComps.end()) {
                    parent[curAnchor] = v;
                    break;
                }
            }
            assert(parent[curAnchor] != none);
        };

        if (sequence.size() == 2) {
            // Happens when the current component is the root component in the topological
            // order. In this case, the root plays the anchor's role.
            if (curAnchor == none) {
                const node v = (sequence[0] == root) ? sequence[1] : sequence[0];
                assert(sequence.front() == root || sequence.back() == root);
                assert(v != root);
                parent[v] = root;
            } else {
                const node v = (sequence[0] == curAnchor) ? sequence[1] : sequence[0];
                assert(v != curAnchor);
                assert(v != root);
                parent[v] = curAnchor;
                if (parent[curAnchor] == none)
                    updateParentOfAnchor();
            }
#ifdef NETWORKIT_SANITY_CHECKS
            checkTwoNodesSequence(sequence);
#endif
            continue;
        }

        // We start building the spanning tree from the first node of the
        // sequence.
        // Root of the spanning tree
        status[sequence[0]] = NodeStatus::IN_SPANNING_TREE;
        const node curAnchorParent = (curAnchor != none) ? parent[curAnchor] : none;

        // All the remaining nodes in the components need to be visited.
        std::for_each(sequence.begin() + 1, sequence.end(),
                      [&status](node u) { status[u] = NodeStatus::NOT_VISITED; });

        count nodesInSpanningTree = 1;
        // Iterate over the remaining nodes to create the spanning tree.
        for (auto it = sequence.begin() + 1; it != sequence.end(); ++it) {
            const node walkStart = *it;
            if (status[walkStart] == NodeStatus::IN_SPANNING_TREE)
                // Node already added to the spanning tree
                continue;

            node currentNode = walkStart;

            // Start a new random walk from the current node.
            do {
                // Get a random neighbor within the component
                node randomNeighbor;
                do {
                    randomNeighbor = G.getIthNeighbor(currentNode, degDist[currentNode](generator));
                } while (status[randomNeighbor] == NodeStatus::NOT_IN_COMPONENT);

                assert(randomNeighbor != none);
                parent[currentNode] = randomNeighbor;
                currentNode = randomNeighbor;

            } while (status[currentNode] != NodeStatus::IN_SPANNING_TREE);

            // Last node encountered in the random walk (it is in the spanning tree);
            const node walkEnd = currentNode;
            assert(status[walkEnd] == NodeStatus::IN_SPANNING_TREE);
            // Add the random walk to the spanning tree; eventually, reverse the path if the
            // anchor/root is encountered
            for (currentNode = walkStart; currentNode != walkEnd;
                 currentNode = parent[currentNode]) {

                status[currentNode] = NodeStatus::IN_SPANNING_TREE;
                ++nodesInSpanningTree;
                if (currentNode == curAnchor || currentNode == root) {

                    // Anchor of current component in the walk, we have to reverse the
                    // parent pointers
                    node next = parent[currentNode];
                    node nextParent = none;
                    node tmp = currentNode;
                    do {
                        status[next] = NodeStatus::IN_SPANNING_TREE;
                        nextParent = parent[next];
                        parent[next] = currentNode;
                        currentNode = next;
                        next = nextParent;
                    } while (next != none);

                    if (tmp == root)
                        parent[root] = none;
                    else if (parent[curAnchor] == none)
                        // Not none if articulation node visited by the parent component before
                        updateParentOfAnchor();
                    else
                        parent[curAnchor] = curAnchorParent;

                    break;
                }
            }

            if (nodesInSpanningTree == sequence.size())
                break;
        }

        for (node u : sequence)
            status[u] = NodeStatus::NOT_IN_COMPONENT;
    }

    count visitedNodes = 0;
    for (node u : G.nodeRange()) {
        while (status[u] == NodeStatus::NOT_IN_COMPONENT) {
            status[u] = NodeStatus::NOT_VISITED;
            ++visitedNodes;
            node parentU = parent[u];
            if (parent[u] != none) {
                assert(siblingPtr[u] == none);
                if (childPtr[parentU] != none)
                    siblingPtr[u] = childPtr[parentU];
                childPtr[parentU] = u;
                u = parentU;
            } else
                break;
        }
        if (visitedNodes == G.numberOfNodes())
            break;
    }

#ifdef NETWORKIT_SANITY_CHECKS
    checkUST();
#endif
}

void ApproxElectricalCloseness::dfsUST() {
    auto &tVisit = tVisitGlobal[omp_get_thread_num()];
    auto &tFinish = tFinishGlobal[omp_get_thread_num()];
    const auto &childPtr = ustChildPtrGlobal[omp_get_thread_num()];
    const auto &siblingPtr = ustSiblingPtrGlobal[omp_get_thread_num()];

    std::stack<std::pair<node, node>> stack;
    stack.push({root, childPtr[root]});

    count timestamp = 0;
    do {
        // v is a child of u that has not been visited yet.
        const node u = stack.top().first;
        const node v = stack.top().second;

        if (v == none) {
            stack.pop();
            tFinish[u] = ++timestamp;
        } else {
            stack.top().second = siblingPtr[v];
            tVisit[v] = ++timestamp;
            stack.push({v, childPtr[v]});
            assert(parentGlobal[omp_get_thread_num()][v] == u);
        }
    } while (!stack.empty());
}

void ApproxElectricalCloseness::aggregateUST() {
#ifdef NETWORKIT_SANITY_CHECKS
    checkTimeStamps();
#endif

    auto &approxEffResistance = approxEffResistanceGlobal[omp_get_thread_num()];
    const auto &tVisit = tVisitGlobal[omp_get_thread_num()];
    const auto &tFinish = tFinishGlobal[omp_get_thread_num()];
    const auto &parent = parentGlobal[omp_get_thread_num()];

    // Doing aggregation
    G.forNodes([&](const node u) {
        node p = bfsParent[u];
        node c = u;

        auto goUp = [&]() -> void {
            c = p;
            p = bfsParent[p];
        };

        while (p != none) {
            // Edge in BSTree: e1 -> e2
            node e1 = p, e2 = c;
            bool reverse = false;
            if (e1 != parent[e2]) {
                if (e2 != parent[e1]) {
                    goUp();
                    continue;
                }
                std::swap(e1, e2);
                reverse = true;
            }

            if (tVisit[u] >= tVisit[e2] && tFinish[u] <= tFinish[e2])
                approxEffResistance[u] += reverse ? -1 : 1;

            goUp();
        }
    });
}

void ApproxElectricalCloseness::run() {
    // Preprocessing
    computeNodeSequence();
    computeBFSTree();

    const count numberOfUSTs = computeNumberOfUSTs();
    Vector sol(G.numberOfNodes());

#pragma omp parallel
    {
        // Thread 0 solves the linear system
        if (omp_get_thread_num() == 0) {
            const count n = G.numberOfNodes();

            const auto L = CSRMatrix::laplacianMatrix(G);
            Diameter diamAlgo(G, estimatedRange, 0);
            diamAlgo.run();
            // Getting diameter upper bound
            const auto diam = static_cast<double>(diamAlgo.getDiameter().second);
            const double tol =
                epsilon * kappa
                / (std::sqrt(static_cast<double>((n * G.numberOfEdges())) * std::log(n)) * diam
                   * 3.);
            ConjugateGradient<CSRMatrix, DiagonalPreconditioner> cg(tol);
            cg.setupConnected(L);

            Vector rhs(n);
            G.forNodes([&](node u) { rhs[u] = -1.0 / static_cast<double>(n); });
            rhs[root] += 1.;
            cg.solve(rhs, sol);

            // ensure column sum of entries is 0.
            double columnSum = 0.;
            G.forNodes([&](node u) { columnSum += sol[u]; });
            sol -= columnSum / static_cast<double>(n);
        }

        // All threads sample and aggregate USTs in parallel
#pragma omp for
        for (omp_index i = 0; i < numberOfUSTs; ++i) {
            sampleUST();
            dfsUST();
            aggregateUST();
        }
    }

    // Aggregating thread-local results
    G.parallelForNodes([&](const node u) {
        // Accumulate all results on thread 0 vector
        for (count i = 1; i < approxEffResistanceGlobal.size(); ++i)
            approxEffResistanceGlobal[0][u] += approxEffResistanceGlobal[i][u];
        diagonal[u] = static_cast<double>(approxEffResistanceGlobal[0][u])
                      / static_cast<double>(numberOfUSTs);
    });

    G.parallelForNodes([&](node u) { diagonal[u] = diagonal[u] - sol[root] + 2. * sol[u]; });
    diagonal[root] = sol[root];
    const double trace = std::accumulate(diagonal.begin(), diagonal.end(), 0.);
    const auto n = static_cast<double>(G.numberOfNodes());

    G.parallelForNodes([&](node u) { scoreData[u] = (n - 1.) / (n * diagonal[u] + trace); });

    hasRun = true;
}

std::vector<double> ApproxElectricalCloseness::computeExactDiagonal(double tol) const {
    Lamg<CSRMatrix> lamg(tol);
    lamg.setupConnected(CSRMatrix::laplacianMatrix(G));

    const count n = G.numberOfNodes();
    const count maxThreads = static_cast<count>(omp_get_max_threads());

    // Solution vectors: one per thread
    std::vector<Vector> solutions(maxThreads, Vector(n));

    // Right hand side vectors: one per thread
    std::vector<Vector> rhss(maxThreads, Vector(n));

    std::vector<double> diag(n);

    const count iters = (n % maxThreads == 0) ? n / maxThreads : n / maxThreads + 1;
    for (count i = 0; i < iters; ++i) {
        // Index of the next vertex to process
        const index base = i * maxThreads;

#pragma omp parallel
        {
            // Each thread solves a linear system from `base` to `base + #threads - 1`
            const index thread = omp_get_thread_num();
            const node v = base + thread;
            if (v < n) {
                // Reset solution and rhs vector of the current thread
                solutions[thread].fill(0.0);

                // Set up system to compute the diagonal entry L^+[v, v]
                rhss[thread].fill(-1. / static_cast<double>(n));
                rhss[thread][v] += 1.;
            }
        }

        if (base + maxThreads >= n) {
            // Last iteration: some threads cannot be used.
            // Resize rhss and solutions to the number of vertices left to be processed.
            rhss.resize(n - base);
            solutions.resize(rhss.size());
        }

        lamg.parallelSolve(rhss, solutions);

        // Store the results
        for (index idx = 0; idx < solutions.size(); ++idx) {
            const node v = base + idx;
            if (v < n)
                diag[v] = solutions[idx][v];
            else
                break;
        }
    }

    return diag;
}

#ifdef NETWORKIT_SANITY_CHECKS
/*
 * Methods for sanity check.
 */
void ApproxElectricalCloseness::checkUST() const {
    std::vector<bool> visitedNodes(G.upperNodeIdBound());
    const auto &parent = parentGlobal[omp_get_thread_num()];
    // To debug
    G.forNodes([&](node u) {
        if (u == root) {
            assert(parent[u] == none);
        } else {
            assert(parent[u] != none);
            std::fill(visitedNodes.begin(), visitedNodes.end(), 0);
            visitedNodes[u] = 1;
            do {
                u = parent[u];
                assert(!visitedNodes[u]);
                visitedNodes[u] = 1;
            } while (u != root);
        }
    });
}

void ApproxElectricalCloseness::checkBFSTree() const {
    G.forNodes([&](node u) {
        if (u == root) {
            assert(bfsParent[u] == none);
        } else {
            std::vector<bool> visited(G.upperNodeIdBound());
            visited[u] = true;
            do {
                u = bfsParent[u];
                assert(!visited[u]);
                assert(!visited[u]);
            } while (u != root);
        }
    });
}

void ApproxElectricalCloseness::checkTwoNodesSequence(const std::vector<node> &sequence) const {
    const std::vector<node> &parent = parentGlobal[omp_get_thread_num()];
    for (node u : sequence) {
        if (u == root) {
            assert(parent[u] == none);
        } else {
            std::vector<bool> visited(G.upperNodeIdBound());
            visited[u] = true;
            do {
                u = parent[u];
                assert(!visited[u]);
                visited[u] = true;
            } while (u != root);
        }
    }
}

void ApproxElectricalCloseness::checkTimeStamps() const {
    const auto &tVisit = tVisitGlobal[omp_get_thread_num()];
    const auto &tFinish = tFinishGlobal[omp_get_thread_num()];
    G.forNodes([&](const node u) {
        assert(tVisit[u] < tFinish[u]);
        if (u == root)
            assert(tVisit[u] == 0);
        else
            assert(tVisit[u] > 0);
        assert(tVisit[u] < 2 * G.numberOfNodes());
    });
}

#endif // NETWORKIT_SANITY_CHECKS

} // namespace NetworKit
