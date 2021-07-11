// no-networkit-format
/*
 * ParallelAgglomerativeClusterer.cpp
 *
 *  Created on: 30.10.2012
 *      Author: Christian Staudt,
 *              Henning Meyerhenke
 */

#include <networkit/coarsening/ClusteringProjector.hpp>
#include <networkit/coarsening/MatchingCoarsening.hpp>
#include <networkit/community/ParallelAgglomerativeClusterer.hpp>
#include <networkit/matching/PathGrowingMatcher.hpp>
#include <networkit/scoring/ModularityScoring.hpp>

namespace NetworKit {

ParallelAgglomerativeClusterer::ParallelAgglomerativeClusterer(const Graph& G) : CommunityDetectionAlgorithm(G) {}

void ParallelAgglomerativeClusterer::run() {

    count MIN_NUM_COMMUNITIES = 2;
    double REL_REPEAT_THRSH = 5e-3; ///< threshold for minimum number of matching edges relative to number of vertices to proceed agglomeration

    // copy graph because we make changes due to merges
    Graph Gcopy(G->numberOfNodes(), true); // make weighted copy
    G->forEdges([&](node u, node v, edgeweight w){
        Gcopy.addEdge(u, v, w);
    });

    std::vector<std::vector<node> > mapHierarchy;

    bool repeat = true;
    do {
        // prepare attributes for scoring
        // FIXME: update to new edge attribute system
        //int attrId = Gcopy.addEdgeAttribute_double(0.0);
        int attrId = 0;

        // perform scoring
        TRACE("before scoring graph of size " , Gcopy.numberOfNodes());
        ModularityScoring<double> modScoring(Gcopy);
        modScoring.scoreEdges(attrId);

        // FIXME: so far only sequential
        // compute matching
        PathGrowingMatcher parMatcher(Gcopy);
        parMatcher.run();
        Matching M = parMatcher.getMatching();

        // contract graph according to matching, TODO: (and star-like structures)
        MatchingCoarsening matchingContracter(Gcopy, M);
        matchingContracter.run();
        Graph Gcombined = matchingContracter.getCoarseGraph();

        // determine if it makes sense to proceed
        count n = Gcopy.numberOfNodes();
        count cn = Gcombined.numberOfNodes();
        count diff = n - cn;
        repeat = ((diff > 0) &&
                (cn >= MIN_NUM_COMMUNITIES) &&
                ((double) diff / (double) n > REL_REPEAT_THRSH)
                ); // TODO: last condition: no community becomes too big

        // prepare next iteration if there is one
        if (repeat) {
            Gcopy = Gcombined;
            mapHierarchy.push_back(matchingContracter.getFineToCoarseNodeMapping());
            TRACE("Repeat agglomeration with graph of size " , Gcopy.numberOfNodes());
        }
    } while (repeat);

    // vertices of coarsest graph are the clusters
    count cn = Gcopy.numberOfNodes();
    Partition zetaCoarse(cn);
    zetaCoarse.allToSingletons();

    // project clustering back to finest graph
    ClusteringProjector projector;
    Partition zeta = projector.projectBackToFinest(zetaCoarse, mapHierarchy,
           *G);
    result = std::move(zeta);
    hasRun = true;
}

std::string ParallelAgglomerativeClusterer::toString() const {
    std::stringstream strm;
    strm << "ParallelAgglomerativeClusterer";
    return strm.str();
}

} /* namespace NetworKit */
