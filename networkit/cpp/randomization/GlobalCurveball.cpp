/*
 * GlobalCurveball.cpp
 *
 *  Created on: 26.05.2018
 *      Author: Manuel Penschuck <networkit@manuel.jetzt>
 */

#include "GlobalCurveballImpl.hpp"
#include <networkit/auxiliary/Random.hpp>
#include <networkit/randomization/DegreePreservingShuffle.hpp>
#include <networkit/randomization/GlobalCurveball.hpp>

namespace NetworKit {

GlobalCurveball::GlobalCurveball(const Graph &G, count number_of_global_trades, bool allowSelfLoops,
                                 bool degreePreservingShufflePreprocessing)
    : impl(new CurveballDetails::GlobalCurveballImpl{G, allowSelfLoops}),
      numGlobalTrades(number_of_global_trades),
      degreePreservingShuffle(degreePreservingShufflePreprocessing) {
    if (allowSelfLoops && !G.isDirected()) {
        throw std::runtime_error("Self loops are only supported for directed graphs");
    }

    if (!allowSelfLoops && G.numberOfSelfLoops()) {
        throw std::runtime_error("Self loops are forbidden but input graph contains some");
    }

    if (G.isWeighted()) {
        throw std::runtime_error("GlobalCurveball supports only unweighted graphs");
    }
}

// We have to define a "default" destructor here, since the definition of
// CurveballDetails::GlobalCurveballImpl is not known in the header file
GlobalCurveball::~GlobalCurveball() = default;

void GlobalCurveball::run() {
    if (hasRun)
        throw std::runtime_error("Call run method only once");

    auto &prng = Aux::Random::getURNG();

    const std::vector<node> *permutation = nullptr;
    DegreePreservingShuffle dps(impl->getInputGraph());

    if (degreePreservingShuffle) {
        dps.run();
        permutation = &dps.getPermutation();
    }

    CurveballDetails::GlobalTradeSequence<CurveballDetails::FixedLinearCongruentialMap<node>> hash{
        impl->getInputGraph().numberOfNodes(), numGlobalTrades, prng};

    if (impl->getInputGraph().isDirected())
        impl->run<true>(hash, permutation);
    else
        impl->run<false>(hash, permutation);

    hasRun = true;
}

Graph GlobalCurveball::getGraph() {
    assureFinished();
    return impl->getGraph();
}

std::string GlobalCurveball::toString() const {
    return "GlobalCurveball";
}

} // namespace NetworKit
