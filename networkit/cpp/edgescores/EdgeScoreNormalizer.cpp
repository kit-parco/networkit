// no-networkit-format
#include <networkit/edgescores/EdgeScoreNormalizer.hpp>

namespace NetworKit {

    template<typename A>
    EdgeScoreNormalizer<A>::EdgeScoreNormalizer(const Graph &G, const std::vector<A> &score, bool invert, double lower, double upper) :
        EdgeScore<double>(G), input(&score), invert(invert), lower(lower), upper(upper) {}

    template<typename A>
    void EdgeScoreNormalizer<A>::run() {
        A minValue = std::numeric_limits< A >::max();
        A maxValue = std::numeric_limits< A >::lowest();

        G->forEdges([&](node, node, edgeid eid) {
            if ((*input)[eid] < minValue) {
                minValue = (*input)[eid];
            }

            if ((*input)[eid] > maxValue) {
                maxValue = (*input)[eid];
            }
        });

        double factor = (upper - lower) / (maxValue - minValue), offset = lower - minValue * factor;

        if (invert) {
            factor *= -1.0;
            offset = upper - minValue * factor;
        }

        scoreData.resize(G->upperEdgeIdBound(), std::numeric_limits<double>::quiet_NaN());

        G->parallelForEdges([&](node, node, edgeid eid) {
            scoreData[eid] = factor * (*input)[eid] + offset;
        });

        hasRun = true;
    }

    template<typename A>
    double EdgeScoreNormalizer<A>::score(node, node) {
        throw std::runtime_error("Not implemented: Use scores() instead.");
    }

    template<typename A>
    double EdgeScoreNormalizer<A>::score(edgeid) {
        throw std::runtime_error("Not implemented: Use scores() instead.");
    }


    template class EdgeScoreNormalizer<double>;
    template class EdgeScoreNormalizer<count>;

} /* namespace NetworKit */
