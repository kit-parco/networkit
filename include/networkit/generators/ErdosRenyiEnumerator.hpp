// no-networkit-format
/*
 * ErdosRenyiGenerator.hpp
 *
 *  Created on: 07.08.2018
 *      Author: Manuel Penschuck (networkit@manuel.jetzt)
 */

#ifndef NETWORKIT_GENERATORS_ERDOS_RENYI_ENUMERATOR_HPP_
#define NETWORKIT_GENERATORS_ERDOS_RENYI_ENUMERATOR_HPP_

#include <omp.h>

#include <atomic>
#include <cassert>
#include <cmath>
#include <random>

#include <networkit/Globals.hpp>
#include <networkit/auxiliary/Random.hpp>
#include <networkit/auxiliary/SignalHandling.hpp>

namespace NetworKit {

/**
 * Generates a stream of edges of a G(n,p) graph. The edges are not
 * written to memory, but emitted only via usual @a forEdges semantics
 * to a callback.
 *
 * Use @ref ErdosRenyiGenerator as a wrapper to output a graph.
 *
 * The enumerator implements both floating point and fixed point arithmetic,
 * to compute the edges which can be selected via the template parameter
 * @a UseFixedPoint. It defaults to true, as this algorithm is typically
 * 2 to 3 times faster. In theory, floating point arithmetic allows for
 * longer runs of consecutive edges not selected. As those events virtually
 * never occur, there are no measurable implications of using the faster
 * variant.
 */
template <bool UseFixedPoint = true>
class ErdosRenyiEnumerator final {
    //! this type is used only internally for fixed-point arithmetic
    using integral_t = unsigned long long;

public:
    /**
     * Generates a G(n, p) graph for n > 1 and 0 < p < 1.
     *
     * For an @b directed graph, the resulting edge stream is equivalent
     * to throwing a coin for each node pair (u, v) and accepting it
     * with probability @a prob; i.e. the stream may include self-loops.
     * Hence the expected number of edges is n*n*prob.
     *
     * For an @b undirected graph, all node pairs (u, v) with 1 < v < u < n
     * are considered. Hence the expected number of edges is n*(n-1)/2*prob
     *
     * @param n         Number of nodes to generate
     * @param prob      Probability that an edge exists
     * @param directed  Selects an directed graph
     */
    ErdosRenyiEnumerator(node n, double prob, bool directed) :
        n{n},
        prob{prob},
        directed{directed}
    {
        assert(n > 0);
    }

    /**
     * Generates an Erdos-Renyi-Graph as specified in the constructor on the fly.
     * The stream is generated in parallel on as many core as an OpenMP parallel
     * section is allotted to. For each edge the callback @a handle is invoked.
     * Two signatures are supported for callback:
     *
     *   (unsigned tid, node u, node v)
     *   (node u, node v)
     *
     * where tid is the current thread id as returned by omp_get_thread_num() and
     * u and v are the edge's nodes. In case of an undirected graph u > v.
     *
     * It is guaranteed that no two threads emit edges for the same u.
     *
     * It can be expected that all threads emit a similar number of edges.
     *
     * Returns number of edges produced.
     */
    template<typename Handle>
    count forEdgesParallel(Handle handle) {
        std::atomic<count> numEdges {0};
        if (directed) {
            #pragma omp parallel
            {
                const unsigned threads = omp_get_num_threads();
                const unsigned tid = omp_get_thread_num();

                const node chunk_size = (n + threads - 1) / threads;
                const node first_node = std::min<node>(n, tid * chunk_size);
                const node last_node  = std::min<node>(n, (tid+1) * chunk_size);

                const auto localNumEdges = enumerate<true>(handle, tid, prob, first_node, last_node);
                numEdges.fetch_add(localNumEdges, std::memory_order_relaxed);
            }
        } else {
            #pragma omp parallel
            {
                const unsigned threads = omp_get_num_threads();
                const unsigned tid = omp_get_thread_num();

                // cells in adj matrix per thread
                const node chunk_size = (n * (n-1) / 2 + threads - 1) / threads;

                double first_node;
                double last_node = 0.0;

                for(unsigned i = 0; i <= tid; i++) {
                    first_node = last_node;
                    node upper_node = std::ceil(std::sqrt(
                        0.25 + first_node * first_node + first_node + 2*chunk_size));
                    last_node = std::min<node>(n, upper_node);
                }

                if (tid + 1 == threads) last_node = n;

                if (first_node < last_node) {
                    const auto localNumEdges = enumerate<false>(handle, tid, prob, first_node, last_node);
                    numEdges.fetch_add(localNumEdges, std::memory_order_relaxed);
                }
            }
        }

        return numEdges.load();
    }

    /**
     * Similarly to @ref forEdgesParallel but computed on one thread only.
     * If the callback accepts three arguments tid is always 0.
     */
    template<typename Handle>
    count forEdges(Handle handle) {
        if (directed) {
            return enumerate<true>(handle, 0, prob, 0, n);
        } else {
            return enumerate<false>(handle, 0, prob, 0, n);
        }
    }

    /**
     * Returns the expected number of edges to be generated.
     */
    count expectedNumberOfEdges() const {
        return prob * n * (directed ? n : (n-1) * 0.5);
    }


private:
    const node n; //< number of nodes
    const double prob; //< probability p
    const bool directed; //< true if a directed graph should be generated

// In the undirected case we only traverse the lower triangle (excluding the
// diagonal) of the adjacency matrix
    template <bool Directed, typename Handle>
    count enumerate(Handle handle, unsigned tid, double prob, const node node_begin, const node node_end) const {
        if (prob > 0.9) {
            // for p > 0.5 we invert the generator and draw the edges NOT in the graph.
            // While this does not change the asymptotical work, it decrease the
            // random bits drawn from prng.

            node cur_u = node_begin;
            node cur_v = 0;
            count num_edges = 0;

            if (!Directed && cur_u == 0)
                cur_u = 1; // all edges need to be of form u > v!

            auto complement_graph = [&] (unsigned tid, node u, node v) {
                while(!(cur_u == u && cur_v == v)) {
                    callHandle(handle, tid, cur_u, cur_v);
                    num_edges++;

                    if (++cur_v == (Directed ? n : cur_u)) {
                        cur_v = 0;
                        cur_u++;
                    }
                }
            };

            enumerate_<Directed>(
                complement_graph, tid, 1.0 - prob, node_begin, node_end);
            complement_graph(tid, node_end, 0);

            return num_edges;
        }

        return enumerate_<Directed>(handle, tid, prob, node_begin, node_end);
    }



    template <bool Directed, typename Handle>
    count enumerate_(Handle handle, unsigned tid, double prob, const node node_begin, const node node_end) const {
        Aux::SignalHandler handler;

        if (prob < std::pow(n, -3.0)) {
            // Even with a union bound over all candidates, WHP no edge will be emitted
            return 0;
        }


        const double inv_log2_cp = 1.0 / std::log2(1.0 - prob);

        // random source
        auto& prng = Aux::Random::getURNG(); // this is thread local
        auto distr = get_distribution<UseFixedPoint>();

        count curr = node_begin;
        if (!Directed && !curr) curr = 1;
        node next = -1;

        node max_skip = 0;
        count numEdges = 0;

        while (curr < node_end) {
            handler.assureRunning();
            // compute new step length
            auto skip = skip_distance(distr(prng), inv_log2_cp);
            next += skip;
            if (skip > max_skip) max_skip = skip;

            // check if at end of row (assuming an average degree of 1,
            // its typically faster to repeatedly subtract than to
            // divide; it is in any case easier easier ...)
            while (next >= (Directed ? n : curr)) {
                // adapt to next row
                next = next - (Directed ? n : curr);
                curr++;
            }

            // insert edge
            if (curr < node_end) {
                numEdges++;
                callHandle(handle, tid, curr, next);
            }
        }

        return numEdges;
    }

// Optimized version of the computation of the skip distance as
// proposed Batagelj and Brandes. It basically converts a uniform
// variate to a geometric random variable.
    count skip_distance(integral_t random_prob, double inv_log2_cp) const {
        /*
         * The original idea is to compute
         * 1 + floor(log(1 - x) / log(1 - prob)) where x is
         * a uniform real from (0, 1).
         *
         * We now precompute inv_log_cp := 1.0 / log(1 - prob) once
         * to avoid recomputing the second log and to replace a
         * division by a multiplication. Hence we compute
         * 1 + floor(log(1 - x) * inv_log_cp).
         *
         * Observe that P[x = k] = P[x = 1-k] and hence we avoid the subtraction:
         * 1 + floor(log(x) * inv_log_cp).
         *
         * Then, typically log() is slower than log2(). On my
         * machines its a factor of roughly 1.5. Thus we replace
         * log by log2 in both equations:
         *  inv_log2_cp := 1.0 / log2(1 - prob)
         *  1 + floor(log2(x) * inv_log2_cp)
         *
         * Now we optimise the generation of the random number.
         * uniform_real_distribution is roughly 3 times slower than
         * uniform_int_distribution. Hence let's switch to fix-point arithmetic.
         * Let X a real drawn uniformly from (0, 2**64), i.e. we model
         * X = (2**64) * x:
         *
         *    1 + floor(log2(x) * inv_log2_cp)
         *  = 1 + floor(log2(X/2**64) * inv_log2_cp)
         *  = 1 + floor((log2(X) - 64) * inv_log2_cp).
         */
        return 1 + static_cast<count>(
            floor((log2(random_prob) - 8*sizeof(integral_t)) * inv_log2_cp)
        );
    }

    count skip_distance(double random_prob, double inv_log2_cp) const {
        return 1 + static_cast<count>(floor((log2(random_prob)) * inv_log2_cp));
    }

// SFINAE to determine and construct the right uniform distribution
    template <bool FixedPoint>
    auto get_distribution() const -> typename std::enable_if<FixedPoint, std::uniform_int_distribution<integral_t>>::type {
        return std::uniform_int_distribution<integral_t>{1, std::numeric_limits<integral_t>::max()};
    }

    template <bool FixedPoint>
    auto get_distribution() const -> typename std::enable_if<!FixedPoint, std::uniform_real_distribution<double>>::type {
        return std::uniform_real_distribution<double>{std::nextafter(0.0, 1.0), std::nextafter(1.0, 0.0)};
    }

// SFINAE to allow using functors with and without thread id as parameter
    template <typename Handle>
    auto callHandle(Handle h, unsigned tid, node u, node v) const -> decltype(h(0u, node{0}, node{0})) {
        return h(tid, u, v);
    }

    template <typename Handle>
    auto callHandle(Handle h, unsigned /*tid*/, node u, node v) const -> decltype(h(node{0}, node{0})) {
        return h(u, v);
    }
};

using ErdosRenyiEnumeratorDefault = ErdosRenyiEnumerator<true>;

}

#endif // NETWORKIT_GENERATORS_ERDOS_RENYI_ENUMERATOR_HPP_
