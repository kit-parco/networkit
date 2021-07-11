/*
 * GlobalCurveballImpl.hpp
 *
 *  Created on: 26.05.2018
 *      Author: Manuel Penschuck <networkit@manuel.jetzt>
 */
#ifndef RANDOMIZATION_GLOBAL_CURVEBALL_IMPL_H_
#define RANDOMIZATION_GLOBAL_CURVEBALL_IMPL_H_

#include <algorithm>
#include <cassert>
#include <type_traits>
#include <utility>
#include <vector>

#include <tlx/algorithm/random_bipartition_shuffle.hpp>
#include <tlx/container/radix_heap.hpp>
#include <tlx/unused.hpp>

#include <networkit/auxiliary/Log.hpp>
#include <networkit/auxiliary/SignalHandling.hpp>
#include <networkit/auxiliary/Timer.hpp>
#include <networkit/base/Algorithm.hpp>
#include <networkit/graph/Graph.hpp>
#include <networkit/graph/GraphBuilder.hpp>
#include <networkit/randomization/GlobalCurveball.hpp>

#include "GlobalTradeSequence.hpp"

namespace NetworKit {
namespace CurveballDetails {

template <typename T1, typename T2>
struct PairFirst {
    T1 operator()(const std::pair<T1, T2> &p) { return p.first; }
};

/**
 * Implementation of EM-GCB ("Parallel and I/O-efficient Randomisation of
 * Massive Networks using Global Curveball Trades", CJ Carstens et al., ESA 2018).
 * This class should not be used directly but rather through the wrapper class
 * GlobalCurveball.
 */
class GlobalCurveballImpl {
    using edgelist_type = std::vector<std::pair<node, node>>;
    using extract_type = PairFirst<node, node>;
    using tfp_queue_type = tlx::RadixHeap<std::pair<node, node>, extract_type, node, 256>;

public:
    GlobalCurveballImpl(const Graph &G, bool allowSelfLoops)
        : inputGraph(G), allowSelfLoops(allowSelfLoops) {
        assert(!allowSelfLoops || inputGraph.isDirected());
    }

    template <bool Directed, typename TradeSequence>
    void run(TradeSequence &trade_sequence, const std::vector<node> *permutation = nullptr) {
        Aux::SignalHandler handler;

        if (hasRun) {
            throw std::runtime_error{"Cannot invoke run several times"};
        }

        Aux::Timer timer;
        timer.start();

        std::vector<node> neighbourhood_of_u, neighbourhood_of_v, disjoint_neighbours,
            common_neighbours;

        auto &urng = Aux::Random::getURNG();

        // copy input graph into queue
        tfp_queue_type current_pq;

        {
            Aux::Timer loadTimer;
            loadTimer.start();

            auto permute = [permutation](node u) -> node {
                return permutation ? (*permutation)[u] : u;
            };

            inputGraph.forNodes([&](node orig_u) {
                const auto u = permute(orig_u);
                const auto hashed_u = trade_sequence.hash(u);
                const auto hint_u = current_pq.get_bucket_key(hashed_u);
                inputGraph.forNeighborsOf(orig_u, [&](node v) {
                    v = permute(v);
                    if (!Directed && u > v)
                        return; // only one message per undirected edge

                    const auto hashed_v = trade_sequence.hash(v);

                    if (Directed || hashed_u < hashed_v) {
                        // u is processed before v
                        current_pq.emplace_in_bucket(hint_u, hashed_u, v);
                    } else {
                        current_pq.emplace(hashed_v, hashed_v, u);
                    }
                });
            });

            loadTimer.stop();
            DEBUG("Loading graph took ", loadTimer.elapsedMilliseconds(), "ms.");
        }

        tfp_queue_type next_pq;

        typename tfp_queue_type::bucket_data_type pq_bucket;
        auto receive_neighbours = [&pq_bucket, &current_pq, this](std::vector<node> &neighbourhood,
                                                                  node x) {
            if (current_pq.empty()) {
                neighbourhood.clear();

            } else {
                current_pq.swap_top_bucket(pq_bucket);
                neighbourhood.resize(pq_bucket.size());
                std::transform(pq_bucket.cbegin(), pq_bucket.cend(), neighbourhood.begin(),
                               [](const std::pair<node, node> &p) { return p.second; });

                assert(inputGraph.degree(x) - neighbourhood.size() <= 1);
                tlx::unused(x, inputGraph);
                pq_bucket.clear();
            }
        };

        for (size_t round = 0; round < trade_sequence.numberOfRounds(); round++) {
            assert(next_pq.empty());
            assert(current_pq.size() == inputGraph.numberOfEdges());
            trade_sequence.switchToRound(round);

            count trade = 0;
            while (!current_pq.empty()) {
                handler.assureRunning();
                trade++;

                // fetch and prepare information received via TFP
                // fetch all messages addressed to next node pair
                const node u = trade_sequence.invert(current_pq.top().first);
                receive_neighbours(neighbourhood_of_u, u);

                node v, hashed_v;

                if (current_pq.empty()) {
                    hashed_v = std::numeric_limits<node>::max();
                    v = u;
                } else {
                    hashed_v = current_pq.top().first;
                    v = trade_sequence.invert(hashed_v);
                }

                receive_neighbours(neighbourhood_of_v, v);

                // filter out edge with trade partner (if it exists, it is
                // contained in u's neighbourhood due to it smaller hash value)
                bool edge_between_uv = false;
                bool edge_between_vu = false;
                if (!Directed || !allowSelfLoops) {
                    auto it = std::find(neighbourhood_of_u.begin(), neighbourhood_of_u.end(), v);
                    if (it != neighbourhood_of_u.end()) {
                        *it = neighbourhood_of_u.back();
                        neighbourhood_of_u.pop_back();
                        edge_between_uv = true;
                    }

                    if (Directed) {
                        auto it =
                            std::find(neighbourhood_of_v.begin(), neighbourhood_of_v.end(), u);
                        if (it != neighbourhood_of_v.end()) {
                            *it = neighbourhood_of_v.back();
                            neighbourhood_of_v.pop_back();
                            edge_between_vu = true;
                        }
                    }
                }

                // compute hashed values of next round
                const auto next_hashed_u = trade_sequence.hashNext(u);
                const auto next_hashed_v = trade_sequence.hashNext(v);

                const auto next_bucket_u = next_pq.get_bucket_key(next_hashed_u);
                const auto next_bucket_v = next_pq.get_bucket_key(next_hashed_v);

                // Split neighborhoods into common and disjoint neighbors
                // Directly forward common neighbours
                disjoint_neighbours.clear();
                common_neighbours.clear();

                const auto num_neighbourhood_of_u = neighbourhood_of_u.size();
                const auto num_neighbourhood_of_v = neighbourhood_of_v.size();

                if (num_neighbourhood_of_u < num_neighbourhood_of_v) {
                    computeCommonDisjointNeighbour(neighbourhood_of_u, neighbourhood_of_v,
                                                   common_neighbours, disjoint_neighbours);
                } else {
                    computeCommonDisjointNeighbour(neighbourhood_of_v, neighbourhood_of_u,
                                                   common_neighbours, disjoint_neighbours);
                }

                // Directly forward neighbours shared by both nodes
                for (const auto neighbour : common_neighbours) {
                    if (Directed) {
                        next_pq.emplace_in_bucket(next_bucket_u, next_hashed_u, neighbour);
                        next_pq.emplace_in_bucket(next_bucket_v, next_hashed_v, neighbour);

                    } else {
                        auto neighbour_hash = trade_sequence.hash(neighbour);
                        if (neighbour_hash > hashed_v) {
                            // neighbour will be processed later this round
                            const auto hint = current_pq.emplace(neighbour_hash, neighbour_hash, u);
                            current_pq.emplace_in_bucket(hint, neighbour_hash, v);

                        } else {
                            // neighbour will be processed in next round,
                            // so send edge into new PQ
                            neighbour_hash = trade_sequence.hashNext(neighbour);

                            const bool u_larger = neighbour_hash < next_hashed_u;
                            const bool v_larger = neighbour_hash < next_hashed_v;

                            if (u_larger && v_larger) {
                                const auto hint =
                                    next_pq.emplace(neighbour_hash, neighbour_hash, u);
                                next_pq.emplace_in_bucket(hint, neighbour_hash, v);

                            } else {
                                if (u_larger) {
                                    next_pq.emplace(neighbour_hash, neighbour_hash, u);
                                } else {
                                    next_pq.emplace_in_bucket(next_bucket_u, next_hashed_u,
                                                              neighbour);
                                }

                                if (v_larger) {
                                    next_pq.emplace(neighbour_hash, neighbour_hash, v);
                                } else {
                                    next_pq.emplace_in_bucket(next_bucket_v, next_hashed_v,
                                                              neighbour);
                                }
                            }
                        }
                    }
                }

                // Shuffle and send disjoint edges
                {
                    auto send_edge = [&](node trade_node, node trade_node_hash,
                                         size_t trade_node_bucket, node neighbour) {
                        if (Directed) {
                            next_pq.emplace_in_bucket(trade_node_bucket, trade_node_hash,
                                                      neighbour);
                            return;
                        }

                        auto neighbour_hash = trade_sequence.hash(neighbour);
                        if (neighbour_hash > hashed_v) {
                            current_pq.emplace(neighbour_hash, neighbour_hash, trade_node);
                        } else {
                            neighbour_hash = trade_sequence.hashNext(neighbour);

                            if (neighbour_hash < trade_node_hash) {
                                next_pq.emplace(neighbour_hash, neighbour_hash, trade_node);
                            } else {
                                next_pq.emplace_in_bucket(trade_node_bucket, trade_node_hash,
                                                          neighbour);
                            }
                        }
                    };

                    const size_t u_setsize = num_neighbourhood_of_u - common_neighbours.size();
                    const size_t v_setsize = num_neighbourhood_of_v - common_neighbours.size();
                    const size_t setsize = u_setsize + v_setsize;
                    assert(u_setsize + v_setsize == disjoint_neighbours.size());

                    tlx::random_bipartition_shuffle(disjoint_neighbours.begin(),
                                                    disjoint_neighbours.end(), u_setsize, urng);

                    size_t i = 0;
                    for (; i < u_setsize; i++) {
                        send_edge(u, next_hashed_u, next_bucket_u, disjoint_neighbours[i]);
                    }

                    for (; i < setsize; i++) {
                        send_edge(v, next_hashed_v, next_bucket_v, disjoint_neighbours[i]);
                    }
                }

                // Do not forget edge between u and v
                if (edge_between_uv) {
                    if (Directed || next_hashed_u < next_hashed_v) {
                        next_pq.emplace_in_bucket(next_bucket_u, next_hashed_u, v);
                    } else {
                        next_pq.emplace_in_bucket(next_bucket_v, next_hashed_v, u);
                    }
                }

                if (Directed && edge_between_vu) {
                    next_pq.emplace_in_bucket(next_bucket_v, next_hashed_v, u);
                }

                assert(current_pq.size() + next_pq.size() == inputGraph.numberOfEdges());
            }

#ifndef NDEBUG
            // After a global trade at most one node may remain with messages
            if (!current_pq.empty()) {
                current_pq.swap_top_bucket(pq_bucket);
                pq_bucket.clear();
                assert(current_pq.empty());
            }
#endif
            current_pq.clear();
            std::swap(current_pq, next_pq);
        }

        hasRun = true;
        prioQueue = std::move(current_pq);

        timer.stop();
        DEBUG("Trading took ", timer.elapsedMilliseconds(), " milliseconds.");
    }

    Graph getGraph() {
        const bool is_directed = inputGraph.isDirected();
        GraphBuilder builder(inputGraph.numberOfNodes(), false, is_directed);

        if (is_directed) {
            for (; !prioQueue.empty(); prioQueue.pop()) {
                const auto top = prioQueue.top();
                assert(allowSelfLoops || top.first != top.second);
                builder.addHalfOutEdge(top.first, top.second);
                builder.addHalfInEdge(top.second, top.first);
            }

        } else {

            for (; !prioQueue.empty(); prioQueue.pop()) {
                const auto top = prioQueue.top();
                builder.addHalfEdge(top.first, top.second);

                if (top.first < top.second)
                    prioQueue.emplace(top.second, top.second, top.first);
            }
        }

        return builder.toGraph(false, true);
    }

    const Graph &getInputGraph() const { return inputGraph; }

protected:
    bool hasRun{false};
    tfp_queue_type prioQueue;

    const Graph &inputGraph;
    const bool allowSelfLoops; ///< Allow self loops (only relevant for directed graphs)

    void computeCommonDisjointNeighbour(std::vector<node> &neighbourhood_of_u,
                                        const std::vector<node> &neighbourhood_of_v,
                                        std::vector<node> &common_neighbours,
                                        std::vector<node> &disjoint_neighbours) const {

        constexpr node BIT = node(1) << (sizeof(node) * 8 - 1);
        constexpr node MASK = ~BIT;

        assert(common_neighbours.empty());
        assert(disjoint_neighbours.empty());

        std::sort(neighbourhood_of_u.begin(), neighbourhood_of_u.end());

        size_t remaining_hits = neighbourhood_of_u.size() / 2;

#ifndef NDEBUG
        const size_t initial_size = neighbourhood_of_u.size() + neighbourhood_of_v.size();
#endif

        for (const auto nv : neighbourhood_of_v) {
            const auto u_it =
                std::lower_bound(neighbourhood_of_u.begin(), neighbourhood_of_u.end(), nv,
                                 [=](const node u, const node v) { return (u & MASK) < v; });

            if (u_it != neighbourhood_of_u.cend() && *u_it == nv) {
                common_neighbours.push_back(nv);
                *u_it |= BIT;

                if (!--remaining_hits) {
                    auto new_end =
                        std::remove_if(neighbourhood_of_u.begin(), neighbourhood_of_u.end(),
                                       [=](const node u) { return u & BIT; });
                    neighbourhood_of_u.resize(std::distance(neighbourhood_of_u.begin(), new_end));
                    remaining_hits = neighbourhood_of_u.size() / 2;
                    if (remaining_hits < 8)
                        remaining_hits = neighbourhood_of_u.size();
                }

            } else {
                disjoint_neighbours.push_back(nv);
            }
        }

        auto new_end = std::remove_if(neighbourhood_of_u.begin(), neighbourhood_of_u.end(),
                                      [=](const node u) { return u & BIT; });
        disjoint_neighbours.insert(disjoint_neighbours.end(), neighbourhood_of_u.begin(), new_end);

        assert(2 * common_neighbours.size() + disjoint_neighbours.size() == initial_size);
#ifndef NDEBUG
        for (auto x : common_neighbours)
            assert(std::find(disjoint_neighbours.cbegin(), disjoint_neighbours.cend(), x)
                   == disjoint_neighbours.cend());
#endif
    }
};

} // namespace CurveballDetails
} // namespace NetworKit

#endif // ! RANDOMIZATION_GLOBAL_CURVEBALL_IMPL_H_
