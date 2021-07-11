// no-networkit-format
/*
 * EdgeScoring.hpp
 *
 *  Created on: 15.10.2012
 *      Author: Christian Staudt
 */

#ifndef NETWORKIT_SCORING_EDGE_SCORING_HPP_
#define NETWORKIT_SCORING_EDGE_SCORING_HPP_

#include <networkit/graph/Graph.hpp>

namespace NetworKit {

/**
 * @ingroup scoring
 * Abstract base class for algorithms associating a score with an edge.
 */
template<typename T>
class EdgeScoring {

protected:
    Graph* G;

public:

    EdgeScoring(Graph& G);

    virtual ~EdgeScoring();

    virtual void scoreEdges(int attrId) = 0;

    virtual T edgeScore(node u, node v) const = 0;
};


template<typename T>
EdgeScoring<T>::EdgeScoring(Graph& G) {
    this->G = &G;
}

template<typename T>
EdgeScoring<T>::~EdgeScoring() {}

} /* namespace NetworKit */
#endif // NETWORKIT_SCORING_EDGE_SCORING_HPP_
