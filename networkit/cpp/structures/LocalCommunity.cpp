#include <networkit/structures/LocalCommunity.hpp>

namespace NetworKit {
template <bool ShellMaintainsExtDeg, bool MaintainBoundary, bool AllowRemoval>
LocalCommunity<ShellMaintainsExtDeg, MaintainBoundary, AllowRemoval>::LocalCommunity(const Graph &G)
    : G(&G), intWeight(0), extWeight(0) {
    if (G.isDirected()) {
        throw std::runtime_error("Directed graphs are not supported");
    }
}

template <bool ShellMaintainsExtDeg, bool MaintainBoundary, bool AllowRemoval>
void LocalCommunity<ShellMaintainsExtDeg, MaintainBoundary, AllowRemoval>::addNode(node u) {
    typename decltype(community)::iterator uIt;
    std::tie(uIt, std::ignore) = community.insert({u, CommunityInfo()});
    shell.erase(u);

    node boundaryNeighbor = none; // if u is in the boundary and has only one neighbor outside of
                                  // the community, store it here.
    std::unordered_map<node, count>::iterator boundaryIt;

    if (MaintainBoundary) {
        boundaryIt = currentBoundary->end();
    }

    G->forNeighborsOf(
        u, [&](node, node v, edgeweight ew) { // insert external neighbors of u into shell
            auto vIt = community.find(v);
            if (vIt != community.end()) {
                if (MaintainBoundary) {
                    auto it = currentBoundary->find(v);
                    assert(it != currentBoundary->end());
                    it->second -= 1;
                    if (it->second == 0) {
                        currentBoundary->erase(it);

                        if (AllowRemoval) {
                            // u was v's only external neighbor
                            *vIt->second.exclusiveOutsideNeighbor = none;

                            // v is now fully internal! -> inform neighbors!
                            G->forNeighborsOf(v, [&](node x) {
                                auto it = community.find(x);
                                assert(it != community.end());

                                it->second.numFullyInternalNeighbors += 1;
                            });
                        }
                    } else if (it->second == 1) {
                        G->forNeighborsOf(v, [&](node x) {
                            auto it = shell.find(x);
                            if (it != shell.end()) {
                                it->second.numExclusiveBoundaryMembers += 1;
                                if (AllowRemoval) {
                                    *uIt->second.exclusiveOutsideNeighbor = it->first;
                                }
                            }
                        });
                    }
                }

                intWeight += ew;
                extWeight -= ew;

                if (AllowRemoval) {
                    uIt->second.intDeg += ew;

                    vIt->second.intDeg += ew;
                    if (ShellMaintainsExtDeg) {
                        vIt->second.extDeg -= ew;
                    }
                }
            } else {
                auto it = shell.find(v);
                if (it == shell.end()) {
                    std::tie(it, std::ignore) = shell.insert({v, ShellInfo()});
                    if (ShellMaintainsExtDeg) {
                        it->second.extDeg.set(G->weightedDegree(v));
                    }
                }

                it->second.intDeg += ew;
                if (ShellMaintainsExtDeg) {
                    it->second.extDeg -= ew;
                }

                extWeight += ew;

                if (AllowRemoval && ShellMaintainsExtDeg) {
                    uIt->second.extDeg += ew;
                }

                if (MaintainBoundary) {
                    if (boundaryIt == currentBoundary->end()) {
                        std::tie(boundaryIt, std::ignore) = currentBoundary->insert({u, 0});
                        boundaryNeighbor = v;
                    }

                    ++boundaryIt->second;
                }

#ifdef NETWORKIT_SANITY_CHECKS
#ifndef NDEBUG
                {
                    auto intExtDeg = calculateIntExtDeg(v);
                    assert(*shell[v].intDeg == intExtDeg.first);
                    if (ShellMaintainsExtDeg) {
                        assert(*shell[v].extDeg == intExtDeg.second);
                    }
                }
#endif
#endif
            }
        });

    if (MaintainBoundary && boundaryIt != currentBoundary->end() && boundaryIt->second == 1) {
        assert(boundaryNeighbor != none);
        shell[boundaryNeighbor].numExclusiveBoundaryMembers += 1;
    }

    if (MaintainBoundary && AllowRemoval && boundaryIt == currentBoundary->end()) {
        // this node is a fully internal node! -> inform neighbors
        G->forNeighborsOf(u, [&](node v) {
            auto it = community.find(v);
            assert(it != community.end());

            it->second.numFullyInternalNeighbors += 1;
        });
    }

#ifdef NETWORKIT_SANITY_CHECKS
#ifndef NDEBUG
    if (AllowRemoval) {
        auto intExtDeg = calculateIntExtDeg(u);
        assert(*community[u].intDeg == intExtDeg.first);
        if (ShellMaintainsExtDeg) {
            assert(*community[u].extDeg == intExtDeg.second);
        }
    }

    assert(!MaintainBoundary || calculateBoundary().size() == currentBoundary->size());
#endif
#endif
}

template <bool ShellMaintainsExtDeg, bool MaintainBoundary, bool AllowRemoval>
void LocalCommunity<ShellMaintainsExtDeg, MaintainBoundary, AllowRemoval>::removeNode(node u) {
    typename decltype(shell)::iterator uIt;
    std::tie(uIt, std::ignore) = shell.insert({u, ShellInfo()});
    community.erase(u);
    bool wasFullyInternal = false;

    if (MaintainBoundary) {
        auto buIt = currentBoundary->find(u);
        if (buIt != currentBoundary->end()) {
            currentBoundary->erase(buIt);
        } else {
            // we are removing a completely internal node!!! (not really a good idea, but okay...)
            wasFullyInternal = true;
        }
    }

    G->forNeighborsOf(
        u, [&](node, node v, edgeweight ew) { // insert external neighbors of u into shell
            auto vIt = community.find(v);
            if (vIt != community.end()) {
                if (MaintainBoundary) {
                    if (wasFullyInternal) {
                        vIt->second.numFullyInternalNeighbors -= 1;
                    }

                    auto it = currentBoundary->find(v);
                    // v was a fully internal node
                    if (it == currentBoundary->end()) {
                        std::tie(it, std::ignore) = currentBoundary->insert({v, 0});
                    }

                    assert(it != currentBoundary->end());

                    it->second += 1;

                    if (it->second == 1) {
                        // v was a fully internal node
                        // therefore, u is now the exclusive outside neighbor of v
                        *vIt->second.exclusiveOutsideNeighbor = u;

                        // inform all neighbors of v that they now have one fully
                        // internal neighbor less
                        G->forNeighborsOf(v, [&](node x) {
                            auto comIt = community.find(x);
                            if (comIt != community.end()) {
                                comIt->second.numFullyInternalNeighbors -= 1;
                            } else {
                                assert(x == u);
                            }
                        });

                        // u has now a neighbor that is only in the boundary
                        // becuase of u
                        uIt->second.numExclusiveBoundaryMembers += 1;
                    } else if (it->second == 2) {
                        *vIt->second.exclusiveOutsideNeighbor = none;

                        uIt->second.numExclusiveBoundaryMembers -= 1;
                    }
                }

                intWeight -= ew;
                extWeight += ew;

                vIt->second.intDeg -= ew;
                uIt->second.intDeg += ew;

                if (ShellMaintainsExtDeg) {
                    vIt->second.extDeg += ew;
                }
            } else {
                auto it = shell.find(v);
                assert(it != shell.end());

                it->second.intDeg -= ew;
                if (ShellMaintainsExtDeg) {
                    it->second.extDeg += ew;
                    uIt->second.extDeg += ew;
                }

                extWeight -= ew;

                if (*it->second.intDeg == 0) {
                    shell.erase(it);
                } else {
#ifdef NETWORKIT_SANITY_CHECKS
#ifndef NDEBUG
                    {
                        auto intExtDeg = calculateIntExtDeg(v);
                        assert(*shell[v].intDeg == intExtDeg.first);
                        if (ShellMaintainsExtDeg) {
                            assert(*shell[v].extDeg == intExtDeg.second);
                        }
                    }
#endif
#endif
                }
            }
        });

#ifdef NETWORKIT_SANITY_CHECKS
#ifndef NDEBUG
    {
        auto intExtDeg = calculateIntExtDeg(u);
        assert(*uIt->second.intDeg == intExtDeg.first);
        if (ShellMaintainsExtDeg) {
            assert(*uIt->second.extDeg == intExtDeg.second);
        }
    }

    assert(!MaintainBoundary || calculateBoundary().size() == currentBoundary->size());
#endif
#endif
}

template <bool ShellMaintainsExtDeg, bool MaintainBoundary, bool AllowRemoval>
bool LocalCommunity<ShellMaintainsExtDeg, MaintainBoundary, AllowRemoval>::contains(node u) const {
    return community.find(u) != community.end();
}

template <bool ShellMaintainsExtDeg, bool MaintainBoundary, bool AllowRemoval>
std::set<node> LocalCommunity<ShellMaintainsExtDeg, MaintainBoundary, AllowRemoval>::toSet() const {
    std::set<node> result;

    for (const auto &it : community) {
        result.insert(it.first);
    }

    return result;
}

template <bool ShellMaintainsExtDeg, bool MaintainBoundary, bool AllowRemoval>
std::unordered_set<node>
LocalCommunity<ShellMaintainsExtDeg, MaintainBoundary, AllowRemoval>::calculateBoundary() {
    std::unordered_set<node> sh;
    for (const auto &it : community) {
        G->forNeighborsOf(it.first, [&](node v) {
            if (!contains(v)) {
                sh.insert(it.first);
            }
        });
    }
    return sh;
}

template <bool ShellMaintainsExtDeg, bool MaintainBoundary, bool AllowRemoval>
std::pair<double, double>
LocalCommunity<ShellMaintainsExtDeg, MaintainBoundary, AllowRemoval>::calculateIntExtDeg(node v) {
    double degInt = 0;
    double degExt = 0;
    G->forNeighborsOf(v, [&](node, node u, edgeweight ew) {
        if (contains(u)) {
            degInt += ew;
        } else {
            degExt += ew;
        }
    });
    return std::make_pair(degInt, degExt);
}

template <bool ShellMaintainsExtDeg, bool MaintainBoundary, bool AllowRemoval>
std::pair<double, double>
LocalCommunity<ShellMaintainsExtDeg, MaintainBoundary, AllowRemoval>::calculateVolumeCut() {
    double internal = 0;
    double external = 0;
    for (const auto &it : community) {
        G->forEdgesOf(it.first, [&](node, node v, edgeweight ew) {
            if (contains(v)) {
                internal += ew;
            } else {
                external += ew;
            }
        });
    }
    internal = internal / 2; // internal edges were counted twice
    return std::make_pair(internal, external);
}

template class LocalCommunity<false, false, false>;
template class LocalCommunity<true, false, false>;
template class LocalCommunity<true, true, false>;
template class LocalCommunity<false, false, true>;
template class LocalCommunity<true, false, true>;
template class LocalCommunity<true, true, true>;

} // namespace NetworKit
