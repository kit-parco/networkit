// no-networkit-format
/*
 * ModularityScoring.hpp
 *
 *  Created on: 15.10.2012
 *      Author: Christian Staudt
 */

#ifndef NETWORKIT_SCORING_MODULARITY_SCORING_HPP_
#define NETWORKIT_SCORING_MODULARITY_SCORING_HPP_

#include <networkit/scoring/EdgeScoring.hpp>
#include <networkit/structures/Partition.hpp>

namespace NetworKit {

// TODO: implement modularity as in Python prototype

/**
 * @ingroup scoring
 */
template<typename T>
class ModularityScoring final : public EdgeScoring<T> {

    double totalEdgeWeight;  //!< total weight of the graph

public:

    /**
     * @param[in]  G  a graph instance
     *
     * Do not modify the graph while using this instance of ModularityScoring.
     */
    ModularityScoring(Graph& G, double gTotalEdgeWeight = 0.0);

    /** Default destructor */
    ~ModularityScoring() override = default;

    void scoreEdges(int attrId) override;

    /**
     * Returns an edge score for an edge (u,v) which expresses the
     * modularity increase which can be gained by merging
     * the clusters of u and v.
     *
     * $$\\Delta mod(c, d) := \\frac{1}{2 \\omega(E)} \\left ( 2 \\omega(E) \\omega(c,d) - \\omega(c) \\omega(d) \\right ) $$
     *
     * @param[in]  u  source node id
     * @param[out]  v  target node id
     *
     */
    T edgeScore(node u, node v) const override;
};

template<typename T>
ModularityScoring<T>::ModularityScoring(Graph& G, double gTotalEdgeWeight) : EdgeScoring<T>(G),
    totalEdgeWeight(gTotalEdgeWeight)
{
    if (gTotalEdgeWeight == 0.0) {
        this->totalEdgeWeight = this->G->totalEdgeWeight();
    }
}

template<typename T>
inline T ModularityScoring<T>::edgeScore(node u, node v) const {
    assert(totalEdgeWeight != 0.0);
    double volume = 2.0 * totalEdgeWeight;
    double nom1 = (this->G->weightedDegree(u) / volume);
    double nom2 = (this->G->weightedDegree(v) / volume);
    double deltaMod = (this->G->weight(u, v) / totalEdgeWeight) -
            (nom1 * nom2);
    return deltaMod;
}

template<typename T>
void ModularityScoring<T>::scoreEdges(int) {
    // TODO: rewrite with new edge attribute system
}

} /* namespace NetworKit */

#endif // NETWORKIT_SCORING_MODULARITY_SCORING_HPP_
