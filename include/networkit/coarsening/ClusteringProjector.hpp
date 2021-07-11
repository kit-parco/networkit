// no-networkit-format
/*
 * ClusteringProjector.hpp
 *
 *  Created on: 07.01.2013
 *      Author: Christian Staudt
 */

#ifndef NETWORKIT_COARSENING_CLUSTERING_PROJECTOR_HPP_
#define NETWORKIT_COARSENING_CLUSTERING_PROJECTOR_HPP_

#include <networkit/graph/Graph.hpp>
#include <networkit/structures/Partition.hpp>

namespace NetworKit {

/**
 * @ingroup coarsening
 */
class ClusteringProjector {

public:

    /**
     * Given
     * 		@param[in]	Gcoarse
     * 		@param[in] 	Gfine
     * 		@param[in]	fineToCoarse
     * 		@param[in]	zetaCoarse	a clustering of the coarse graph
     *
     * 	, project the clustering back to the fine graph to create a clustering of the fine graph.
     * 		@param[out] 			a clustering of the fine graph
     **/
    virtual Partition projectBack(const Graph& Gcoarse, const Graph& Gfine, const std::vector<node>& fineToCoarse, const Partition& zetaCoarse);



    /**
     * Project a clustering \zeta^{i} of the coarse graph G^{i} back to
     * the finest graph G^{0}, using the hierarchy of fine->coarse maps
     */
    virtual Partition projectBackToFinest(const Partition& zetaCoarse, const std::vector<std::vector<node> >& maps, const Graph& Gfinest);


    /**
     * Assuming that the coarse graph resulted from contracting and represents a clustering of the finest graph
     *
     * @param[in]	Gcoarse		coarse graph
     * @param[in]	Gfinest		finest graph
     * @param[in]	maps		hierarchy of maps M^{i->i+1} mapping nodes in finer graph to supernodes in coarser graph
     */
    virtual Partition projectCoarseGraphToFinestClustering(const Graph& Gcoarse, const Graph& Gfinest, const std::vector<std::vector<node> >& maps);

};

} /* namespace NetworKit */
#endif // NETWORKIT_COARSENING_CLUSTERING_PROJECTOR_HPP_
