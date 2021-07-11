/* GCE.cpp
 *
 * Created on: 06.05.2013
 * Author: cls
 */

#include <unordered_map>
#include <utility>

#include <networkit/auxiliary/IncrementalUniformRandomSelector.hpp>
#include <networkit/scd/GCE.hpp>
#include <networkit/structures/LocalCommunity.hpp>

namespace NetworKit {

GCE::GCE(const Graph &g, std::string objective)
    : SelectiveCommunityDetector(g), objective(std::move(objective)) {
    if (g.numberOfSelfLoops()) {
        throw std::runtime_error("Graphs with self-loops are not supported in GCE");
    }
}

std::set<node> GCE::expandSeed(node seed) {
    WARN("GCE::expandSeed is deprecated, use GCE::expandOneCommunity instead");
    return expandOneCommunity(seed);
}

template <bool objectiveIsM>
std::set<node> expandseedInternal(const Graph &g, const std::set<node> &seeds) {
    double currentQ = 0.0; // current community quality

    LocalCommunity<true, !objectiveIsM> com(g);

    for (node s : seeds) {
        com.addNode(s);
#ifdef NETWORKIT_SANITY_CHECKS
        assert(com.contains(s));
#endif
    }

    using shell_info_t = typename decltype(com)::ShellInfo;

    /*
     * objective function M
     * @return quality difference for the move of v to C
     */
    auto deltaM = [&](const shell_info_t &info) {
        double delta = (com.internalEdgeWeight() + (*info.intDeg))
                       / (double)(com.cut() - (*info.intDeg) + (*info.extDeg));
        return delta - currentQ;
    };

    /*
     * objective function L
     * @return quality difference for the move of v to C
     */
    auto deltaL = [&](const shell_info_t &info) {
        // Compute difference in boundary size: for each neighbor where we are the last
        // external neighbor decrease by 1, if v has an external neighbor increase by 1
        int64_t boundaryDiff = info.boundaryChange();

        TRACE("boundary diff: ", boundaryDiff);

        double numerator =
            2.0 * (com.internalEdgeWeight() + (*info.intDeg)) * (com.boundarySize() + boundaryDiff);
        double denominator = (com.size() + 1) * (com.cut() - (*info.intDeg) + (*info.extDeg));
        return (numerator / denominator) - currentQ;
    };

    // select quality objective
    auto deltaQ = [&](const shell_info_t &info) -> double {
        if (objectiveIsM) {
            return deltaM(info);
        } else {
            return deltaL(info);
        }
    };

    if (objectiveIsM) {
        currentQ = com.internalEdgeWeight() / com.cut();
    } else {
        double numerator = 2.0 * com.internalEdgeWeight() * com.boundarySize();
        double denominator = com.size() * com.cut();
        currentQ = numerator / denominator;
    }

    node vMax;
    do {
        // scan shell for node with maximum quality improvement
        double dQMax = 0.0; // maximum quality improvement
        vMax = none;
        Aux::IncrementalUniformRandomSelector selector;

        com.forShellNodes([&](const node v, const shell_info_t &info) {
            // get values for current node
            double dQ = deltaQ(info);

            TRACE("dQ: ", dQ);

            if (dQ > dQMax) {
                vMax = v;
                dQMax = dQ;
                selector.reset();
            } else if (dQ == dQMax && selector.addElement()) {
                vMax = v;
            }
        });

        TRACE("vMax: ", vMax);
        TRACE("dQMax: ", dQMax);
        if (vMax != none) {
            com.addNode(vMax); // add best node to community
            currentQ += dQMax; // update current community quality
        }
    } while (vMax != none);

    return com.toSet();
}

std::set<node> GCE::expandOneCommunity(const std::set<node> &s) {
    if (objective == "M") {
        return expandseedInternal<true>(*g, s);
    } else if (objective == "L") {
        return expandseedInternal<false>(*g, s);
    } else {
        throw std::runtime_error("unknown objective function");
    }
}

} /* namespace NetworKit */
