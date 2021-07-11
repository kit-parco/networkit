// no-networkit-format

#ifndef NETWORKIT_SPARSIFICATION_SCAN_STRUCTURAL_SIMILARITY_SCORE_HPP_
#define NETWORKIT_SPARSIFICATION_SCAN_STRUCTURAL_SIMILARITY_SCORE_HPP_

#include <networkit/edgescores/EdgeScore.hpp>

namespace NetworKit {

class SCANStructuralSimilarityScore final : public EdgeScore<double> {

public:
    SCANStructuralSimilarityScore(const Graph& G, const std::vector<count>& triangles);
    void run() override;
    double score(edgeid eid) override;
    double score(node u, node v) override;

private:
    const std::vector<count>* triangles;

};

} // namespace NetworKit

#endif // NETWORKIT_SPARSIFICATION_SCAN_STRUCTURAL_SIMILARITY_SCORE_HPP_
