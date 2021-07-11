#include <iterator>

#include <networkit/auxiliary/Random.hpp>
#include <networkit/scd/RandomBFS.hpp>

namespace NetworKit {

RandomBFS::RandomBFS(const Graph &g, const Cover &c)
    : SelectiveCommunityDetector(g), c(&c), subsetSizes(c.subsetSizeMap()) {}

std::set<node> RandomBFS::expandOneCommunity(const std::set<node> &s) {
    // allow as many nodes in the community as seeds are given
    count comSize = s.size();

    { // select a random community of s and get its size
        std::set<index> gs(c->subsetsOf(*s.begin()));

        for (auto it = ++s.begin(); it != s.end(); ++it) {
            const std::set<index> additionalComs(c->subsetsOf(*it));
            for (auto it = gs.begin(); it != gs.end();) {
                if (additionalComs.find(*it) == additionalComs.end()) {
                    it = gs.erase(it);
                } else {
                    ++it;
                }
            }
        }

        if (!gs.empty()) {
            index i = Aux::Random::index(gs.size());
            auto it = gs.begin();
            std::advance(it, i);
            comSize = subsetSizes[*it];
        }
    }

    std::set<node> result;

    std::vector<node> currentLevel, nextLevel;

    for (node u : s) {
        currentLevel.push_back(u);
    }

    while (result.size() < comSize && !currentLevel.empty()) {
        nextLevel.clear();

        if (currentLevel.size() + result.size() < comSize) {
            for (node u : currentLevel) {
                result.insert(u);
            }
        } else {
            std::shuffle(currentLevel.begin(), currentLevel.end(), Aux::Random::getURNG());

            for (auto it = currentLevel.begin(); result.size() < comSize; ++it) {
                assert(it != currentLevel.end());
                result.insert(*it);
            }

            break;
        }

        for (node u : currentLevel) {
            g->forNeighborsOf(u, [&](node v) {
                if (result.count(v) == 0) {
                    nextLevel.push_back(v);
                }
            });
        }

        std::sort(nextLevel.begin(), nextLevel.end());
        auto last = std::unique(nextLevel.begin(), nextLevel.end());
        nextLevel.erase(last, nextLevel.end());

        std::swap(currentLevel, nextLevel);
    }

    return result;
}

} // namespace NetworKit
