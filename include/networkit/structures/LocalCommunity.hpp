#ifndef NETWORKIT_STRUCTURES_LOCAL_COMMUNITY_HPP_
#define NETWORKIT_STRUCTURES_LOCAL_COMMUNITY_HPP_

#include <set>
#include <unordered_map>
#include <unordered_set>

#include <networkit/graph/Graph.hpp>

namespace NetworKit {

/**
 * Class for maintaining some measures for a local community while expanding it.
 *
 * Depending on the template parameters, different values are maintained. By
 * default, it stores the shell of the community (i.e., all external neighbors)
 * and the volume and cut of the community. Further, for every node in the shell,
 * the number of internal neighbors is maintained. Optionally, also the number of
 * external neighbors can be maintained for every node in the shell (@a ShellMaintainsExtDeg).
 * Additionally, also the boundary, i.e., all nodes in the community that have
 * a neighbor in the shell are maintained (@a MaintainBoundary).
 */
template <bool ShellMaintainsExtDeg, bool MaintainBoundary = false, bool AllowRemoval = false>
class LocalCommunity {
    static_assert(!MaintainBoundary || ShellMaintainsExtDeg,
                  "boundary cannot be maintained without maintaining the external degree");

public:
    template <typename ValueType, bool used>
    struct OptionalValue {
        ValueType &get() { throw std::runtime_error("Getting value that is missing"); }

        const ValueType &get() const { throw std::runtime_error("Getting value that is missing"); }

        ValueType &operator*() { throw std::runtime_error("Getting value that is missing"); }

        const ValueType &operator*() const {
            throw std::runtime_error("Getting value that is missing");
        }

        ValueType *operator->() { throw std::runtime_error("Getting value that is missing"); }

        const ValueType *operator->() const {
            throw std::runtime_error("Getting value that is missing");
        }

        void set(ValueType) { throw std::runtime_error("Setting value that is missing"); }

        OptionalValue &operator+=(ValueType) {
            throw std::runtime_error("Increasing value that is missing");
        }

        OptionalValue &operator-=(ValueType) {
            throw std::runtime_error("Decreasing value that is missing");
        }

        OptionalValue(const ValueType &){};

        OptionalValue(){};
    };

    template <typename ValueType>
    struct OptionalValue<ValueType, true> {
        ValueType &get() { return value; }

        const ValueType &get() const { return value; }

        ValueType &operator*() { return value; }

        const ValueType &operator*() const { return value; }

        ValueType *operator->() { return &value; }

        const ValueType *operator->() const { return &value; }

        void set(ValueType v) { value = v; }

        template <typename SecondValue,
                  typename std::enable_if<std::is_scalar<SecondValue>::value, int>::type = 0>
        OptionalValue &operator+=(SecondValue v) {
            value += v;
            return *this;
        }

        template <typename SecondValue,
                  typename std::enable_if<std::is_scalar<SecondValue>::value, int>::type = 0>
        OptionalValue &operator-=(SecondValue v) {
            value -= v;
            return *this;
        }

        OptionalValue(const ValueType &v) : value(v){};

        OptionalValue() : value(){};

        ValueType value;
    };

    struct ShellInfo {
        /**
         * Number of neighbors in the community.
         */
        OptionalValue<double, true> intDeg;

        /**
         * Number of neighbors outside the community.
         */
        OptionalValue<double, ShellMaintainsExtDeg> extDeg;

        /**
         * Number of nodes in the boundary whose only outside neighbor is this node.
         */
        OptionalValue<count, MaintainBoundary> numExclusiveBoundaryMembers;

        /**
         * Calculate how much the boundary size would change if the node
         * was added to the community.
         */
        int64_t boundaryChange() const {
            int64_t boundary_diff = -1 * numExclusiveBoundaryMembers.get();

            if (extDeg.get() > 0) {
                boundary_diff += 1;
            }

            return boundary_diff;
        }

        ShellInfo() : intDeg(0), extDeg(0), numExclusiveBoundaryMembers(0) {}
    };

    struct CommunityInfo {
        /**
         * Number of neighbors in the community.
         */
        OptionalValue<double, AllowRemoval> intDeg;

        /**
         * Number of neighbors outside the community.
         */
        OptionalValue<double, ShellMaintainsExtDeg && AllowRemoval> extDeg;

        /**
         * The only neighbor outside the community if there is only one.
         */
        OptionalValue<node, MaintainBoundary && AllowRemoval> exclusiveOutsideNeighbor;

        /**
         * The number of neighbors that have no internal neighbors.
         */
        OptionalValue<count, MaintainBoundary && AllowRemoval> numFullyInternalNeighbors;

        /**
         * Calculate how much the boundary size would change if the node
         * was removed from the community.
         */
        int64_t boundaryChange() const {
            int64_t boundary_diff = numFullyInternalNeighbors.get();

            if (extDeg.get() > 0) {
                boundary_diff -= 1;
            }

            return boundary_diff;
        }

        CommunityInfo()
            : intDeg(0), extDeg(0), exclusiveOutsideNeighbor(none), numFullyInternalNeighbors(0) {}
    };

    /**
     * Initialize an empty community for the given graph.
     *
     * @param G The graph.
     */
    LocalCommunity(const Graph &G);

    /**
     * Add the given node to the community.
     *
     * This needs O(deg(u)) time in expectation. If the boundary
     * is maintained, the complexity can be higher but this is amortized
     * if nodes are only added to the community.
     *
     * @param u The node to add.
     */
    void addNode(node u);

    /**
     * Removes the given node form the community
     *
     * This needs O(deg(u)) time in expectation. If the boundary
     * is maintained, the complexity can be higher but this is amortized
     * if nodes are only added to the community.
     *
     * @param u The node to add.
     */
    void removeNode(node u);

    /**
     * Check if a given node @a u is in the community.
     *
     * @param u The node to check.
     *
     * @return If @a u is in the community.
     */
    bool contains(node u) const;

    /**
     * Return the community as a std::set<node>.
     *
     * @return The community as std::set<node>
     */
    std::set<node> toSet() const;

    /**
     * Execute a callback for every node in the shell.
     *
     * The complexity is linear in the size of the shell.
     * The callback should accept a node and a ShellInfo
     * object as parameters. Note that only the enabled
     * properties of the ShellInfo object may be accessed.
     *
     * @param callback The callback to call.
     */
    template <typename F>
    void forShellNodes(F callback) {
        for (const auto &it : shell) {

#ifdef NETWORKIT_SANITY_CHECKS
#ifndef NDEBUG
            auto intExtDeg = calculateIntExtDeg(it.first);
            assert(*it.second.intDeg == intExtDeg.first);
            if (ShellMaintainsExtDeg) {
                assert(*it.second.extDeg == intExtDeg.second);
            }

            if (MaintainBoundary) {
                int64_t boundary_diff_debug = 0;
                bool v_in_boundary = false;
                G->forNeighborsOf(it.first, [&](node x) {
                    auto it = currentBoundary->find(x);
                    if (it != currentBoundary->end()) {
                        if (it->second == 1) {
                            boundary_diff_debug -= 1;
                        }
                    } else if (!v_in_boundary) {
                        boundary_diff_debug += 1;
                        v_in_boundary = true;
                    }
                });

                assert(it.second.boundaryChange() == boundary_diff_debug);
            }
#endif // NDEBUG
#endif // NETWORKIT_SANITY_CHECKS

            callback(it.first, it.second);
        }
    }

    /**
     * Execute a callback for every node in the community.
     *
     * The complexity is linear in the size of the community.
     * The callback should accept a node and a CommunityInfo
     * object as parameters. Note that only the enabled
     * properties of the CommunityInfo object may be accessed.
     *
     * Note that unless AllowRemoval is set, the CommunityInfo object is empty.
     *
     * @param callback The callback to call.
     */
    template <typename F>
    void forCommunityNodes(F callback) {
        for (auto it = community.cbegin(); it != community.cend();) {
            // advance the iterator so the current element may be deleted
            const auto cit = it++;
#ifdef NETWORKIT_SANITY_CHECKS
#ifndef NDEBUG
            if (AllowRemoval) {
                auto intExtDeg = calculateIntExtDeg(cit->first);
                assert(*cit->second.intDeg == intExtDeg.first);
                if (ShellMaintainsExtDeg) {
                    assert(*cit->second.extDeg == intExtDeg.second);
                }

                if (MaintainBoundary) {
                    int64_t boundary_diff_debug = 0;
                    bool v_in_boundary = false;
                    G->forNeighborsOf(cit->first, [&](node x) {
                        if (community.find(x) == community.end()) {
                            if (!v_in_boundary) {
                                boundary_diff_debug -= 1;
                                v_in_boundary = true;
                            }
                        } else if (currentBoundary->find(x) == currentBoundary->end()) {
                            boundary_diff_debug += 1;
                        }
                    });

                    assert(cit->second.boundaryChange() == boundary_diff_debug);
                }
            }
#endif // NDEBUG
#endif // NETWORKIT_SANITY_CHECKS

            callback(cit->first, cit->second);
        }
    }

    /**
     * Get the size of the community.
     *
     * @return The size of the community.
     */
    size_t size() const { return community.size(); }

    /**
     * Get the internal edge weight.
     * This is the sum of the weights of all internal edges, i.e.,
     * the weight of every internal edge is counted once.
     *
     * @return The internal edge weight of the community.
     */
    double internalEdgeWeight() const { return intWeight; }

    /**
     * Get the cut of the community.
     *
     * @return The cut of the community.
     */
    double cut() const { return extWeight; }

    /**
     * Get the size of the boundary of the community.
     *
     * @return The size of the boundary.
     */
    count boundarySize() const { return currentBoundary->size(); }

private:
    const Graph *G;
    std::unordered_map<node, CommunityInfo> community;
    std::unordered_map<node, ShellInfo> shell;
    double intWeight;
    double extWeight;
    OptionalValue<std::unordered_map<node, count>, MaintainBoundary> currentBoundary;

    // The boundary is defined as all nodes of C that have a neighbor not in C
    std::unordered_set<node> calculateBoundary();

    /**
     * internal and external weighted degree of a node with respect to the community
     */
    std::pair<double, double> calculateIntExtDeg(node v);

    std::pair<double, double> calculateVolumeCut();
};

} // namespace NetworKit

#endif // NETWORKIT_STRUCTURES_LOCAL_COMMUNITY_HPP_
