#include <stdexcept>

#include <networkit/auxiliary/Random.hpp>
#include <networkit/auxiliary/SignalHandling.hpp>
#include <networkit/randomization/DegreePreservingShuffle.hpp>
#include <networkit/randomization/EdgeSwitching.hpp>

namespace NetworKit {

void EdgeSwitchingInPlace::run() {
    auto numberOfSwitches =
        static_cast<count>(std::ceil(graph.numberOfEdges() * numberOfSwitchesPerEdge));

    if (graph.numberOfEdges() < 2 || numberOfSwitches == 0)
        return;

    if (!hasRun) {
        degreeDistribution = std::discrete_distribution<node>(
            graph.numberOfNodes(), //
            0.0, static_cast<double>(graph.numberOfNodes()),
            [&](double x) { return graph.degree(static_cast<node>(x)); });
    }

    auto &urng = Aux::Random::getURNG();
    Aux::SignalHandler handler;
    while (numberOfSwitches--) {
        handler.assureRunning();
        const auto s1 = degreeDistribution(urng);
        const auto s2 = degreeDistribution(urng);

        // we avoid GraphTools::randomNeighbor to avoid the implicit cost of accessing the
        // Aux::Random::getURNG
        const auto i1 = std::uniform_int_distribution<index>{0, graph.degree(s1) - 1}(urng);
        const auto t1 = graph.getIthNeighbor(s1, i1);

        if (s2 == t1 || graph.hasEdge(s2, t1)) // early reject
            continue;

        const auto i2 = std::uniform_int_distribution<index>{0, graph.degree(s2) - 1}(urng);
        const auto t2 = graph.getIthNeighbor(s2, i2);

        if (t1 == t2 || s1 == t2 || graph.hasEdge(s1, t2)) // remaining checks
            continue;

        graph.swapEdge(s1, t1, s2, t2);

        ++numberOfSwapsPerformed;
    }

    hasRun = true;
}

void EdgeSwitchingInPlace::setNumberOfSwitchesPerEdge(double x) {
    if (x < 0)
        throw std::invalid_argument("NumberOfSwitchesPerEdge has to be non negative");

    numberOfSwitchesPerEdge = x;
}

EdgeSwitching::EdgeSwitching(const NetworKit::Graph &G, double numberOfSwitchesPerEdge,
                             bool doPreprocessing)
    : ownedGraph(doPreprocessing ? DegreePreservingShuffle::shuffleGraph(G) : G),
      inPlaceAlgorithm(ownedGraph, numberOfSwitchesPerEdge) {}

} // namespace NetworKit
