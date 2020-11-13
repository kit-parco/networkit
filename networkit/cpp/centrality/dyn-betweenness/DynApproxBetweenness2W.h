/*
 * DynApproxBetweenness2.h
 *
 *  Created on: 16.02.2015
 *      Author: ebergamini
 */

#ifndef DYNAPPROXBETW2W_H_
#define DYNAPPROXBETW2W_H_

#include "Centrality.h"
#include "DynCentrality.h"
#include "../dynamics/GraphEvent.h"
#include "../graph/BFSvisit.h"
#include "../graph/DijkstraVisit.h"


#include <math.h>
#include <algorithm>
#include <memory>
#include <omp.h>

namespace NetworKit {

/**
 * @ingroup graph
 * Interface for dynamic approximated betweenness centrality algorithm.
 */
class DynApproxBetweenness2W: public Centrality, public DynCentrality {

public:
    /**
      * The algorithm approximates the betweenness of all vertices so that the scores are
      * within an additive error @a epsilon with probability at least (1- @a delta).
      * The values are normalized by default.
      *
      * @param	G			the graph
      * @param  storePredecessors   keep track of the lists of predecessors?
      * @param	epsilon		maximum additive error
      * @param	delta		probability that the values are within the error guarantee
     */
    DynApproxBetweenness2W(const Graph& G, double epsilon=0.01, double delta=0.1, bool storePredecessors = true);

    /**
     * Runs the static approximated betweenness centrality algorithm on the initial graph.
     */
    void run() override;

    /**
    * Updates the betweenness centralities after a batch of edge insertions on the graph.
    *
    * @param batch The batch of edge insertions.
    */
    void update(const std::vector<GraphEvent>& batch);

    /**
    * Get number of path samples used for last calculation
    */
    count getNumberOfSamples();

private:
    void sampleNewPaths(count, count);

    bool storePreds = false; //TODO add vector with predecessors
    double epsilon; //!< maximum error
    double delta;
    const double c = 0.5; // universal positive constant - see reference in paper
    count vd;
    count r;
    count r2;
    std::vector<DijkstraVisit> sssp;
    std::vector<DijkstraVisit> sssp2;
    std::vector<std::vector<count>> npaths;
    std::vector<edgeweight> wmin;
    std::vector<count> compSize;

    std::vector<count> vis;
    std::vector<node> u;
    std::vector<node> v;
    std::vector<edgeweight> maxDist;
    std::vector<edgeweight> maxDist2;
    std::vector <std::vector<node>> sampledPaths;

    const edgeweight infDist = std::numeric_limits<edgeweight>::max();

};

} /* namespace NetworKit */

#endif /* DYNAPPROXBETW2W_H_ */
