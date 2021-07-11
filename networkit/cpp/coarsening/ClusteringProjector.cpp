// no-networkit-format
/*
 * ClusteringProjector.cpp
 *
 *  Created on: 07.01.2013
 *      Author: Christian Staudt
 */

#include <networkit/auxiliary/Log.hpp>
#include <networkit/coarsening/ClusteringProjector.hpp>

namespace NetworKit {


Partition ClusteringProjector::projectBack(const Graph&, const Graph& Gfine, const std::vector<node>& fineToCoarse,
        const Partition& zetaCoarse) {

    Partition zetaFine(Gfine.upperNodeIdBound());
    zetaFine.setUpperBound(zetaCoarse.upperBound());
    Gfine.forNodes([&](node v) {
        node sv = fineToCoarse[v];
        index cv = zetaCoarse.subsetOf(sv);
        zetaFine.addToSubset(cv,v);
    });

    return zetaFine;
}

Partition ClusteringProjector::projectBackToFinest(const Partition& zetaCoarse,
        const std::vector<std::vector<node> >& maps, const Graph& Gfinest) {
    if (zetaCoarse.numberOfElements() == Gfinest.numberOfNodes()) {
        return zetaCoarse;
    }

    Partition zetaFine(Gfinest.upperNodeIdBound());
    zetaFine.setUpperBound(zetaCoarse.upperBound()); // upper bound for ids in zetaFine must be set to upper bound of zetaCoarse, or modularity assertions fail

    // store temporarily coarsest supernode here
    std::vector<node> tempMap(Gfinest.upperNodeIdBound());

    // initialize to identity
    Gfinest.parallelForNodes([&](node v){
        tempMap[v] = v;
    });

    // find coarsest supernode for each node
    for (auto iter = maps.begin(); iter != maps.end(); ++iter) {
        Gfinest.parallelForNodes([&](node v){
            tempMap[v] = (* iter)[tempMap[v]];
        });
    }


    // set clusters for fine nodes
    Gfinest.parallelForNodes([&](node v) {
        index sc = zetaCoarse[tempMap[v]];
        zetaFine.addToSubset(sc,v);
    });

    return zetaFine;
}

Partition ClusteringProjector::projectCoarseGraphToFinestClustering(const Graph& Gcoarse, const Graph& Gfinest, const std::vector<std::vector<node> >& maps) {

    Partition zeta(Gfinest.upperNodeIdBound());
    zeta.setUpperBound(Gcoarse.upperNodeIdBound());

    // store temporarily coarsest supernode here
    std::vector<node> super(Gfinest.upperNodeIdBound());
    // initialize to identity
    Gfinest.parallelForNodes([&](node v){
        super[v] = v;
    });


    // find coarsest supernode for each node
    for (auto iter = maps.begin(); iter != maps.end(); ++iter) {
        Gfinest.parallelForNodes([&](node v){
            super[v] = (* iter)[super[v]];
        });
    }

    // assign super node id as cluster id
    Gfinest.parallelForNodes([&](node v) {
        zeta[v] = super[v];
    });

    DEBUG("number of clusters in projected clustering: " , zeta.numberOfSubsets());
    DEBUG("number of nodes in coarse graph: " , Gcoarse.numberOfNodes());
    assert (zeta.numberOfSubsets() == Gcoarse.numberOfNodes());

    return zeta;

}

} /* namespace NetworKit */
