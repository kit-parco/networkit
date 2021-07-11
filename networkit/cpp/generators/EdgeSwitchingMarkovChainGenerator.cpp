
#include <networkit/generators/EdgeSwitchingMarkovChainGenerator.hpp>
#include <networkit/generators/HavelHakimiGenerator.hpp>
#include <networkit/randomization/EdgeSwitching.hpp>

namespace NetworKit {

EdgeSwitchingMarkovChainGenerator::EdgeSwitchingMarkovChainGenerator(
    const std::vector<count> &sequence, bool ignoreIfNotRealizable, count numSwitchesPerEdge)
    : StaticDegreeSequenceGenerator(sequence), ignoreIfNotRealizable(ignoreIfNotRealizable),
      numSwitchesPerEdge(numSwitchesPerEdge) {}

Graph EdgeSwitchingMarkovChainGenerator::generate() {
    Graph result(HavelHakimiGenerator(seq, ignoreIfNotRealizable).generate());

    EdgeSwitchingInPlace edgeSwitching(result, numSwitchesPerEdge);
    edgeSwitching.run();

    return result;
}

} // namespace NetworKit
