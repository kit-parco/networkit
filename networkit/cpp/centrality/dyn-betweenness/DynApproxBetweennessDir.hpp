/*
 * DynApproxBetweennessDir.h
 *
 *  Created on: 31.07.2014
 *      Author: ebergamini
 */

<<<<<<< Updated upstream
<<<<<<< HEAD
#ifndef DYNAPPROXBETW_DIR_H_
#define DYNAPPROXBETW_DIR_H_
=======
#ifndef NETWORKIT_CENTRALITY_DYN_APPROX_BETWEENNESS_DIR_HPP_
#define NETWORKIT_CENTRALITY_DYN_APPROX_BETWEENNESS_DIR_HPP_
>>>>>>> Stashed changes

#include "Centrality.h"
#include "DynCentrality.h"
#include "../dynamics/GraphEvent.h"
#include "../graph/DynSSSP.h"
#include "../graph/BFSvisit.h"

#include <math.h>
=======
/*
 * Note from Charmaine:
 *
 * This class is for insertions+deletions in directed unweighted graphs
 */

#ifndef DYNAPPROXBETW_DIR_H_
#define DYNAPPROXBETW_DIR_H_


#include <networkit/centrality/Centrality.hpp>
#include <networkit/dynamics/GraphEvent.hpp>
#include <networkit/distance/DynSSSP.hpp>

#include <cmath>
>>>>>>> 7ca6f73b0... Fix includes and first steps at compiling files (not compiling yet)
#include <algorithm>
#include <memory>
#include <omp.h>

namespace NetworKit {

/**
 * @ingroup graph
 * Interface for dynamic approximated BetweennessDir centrality algorithm.
 */
<<<<<<< HEAD
class DynApproxBetweennessDir: public Centrality, public DynCentrality {
=======
class DynApproxBetweennessDir: public Centrality, public DynAlgorithm {
>>>>>>> 7ca6f73b0... Fix includes and first steps at compiling files (not compiling yet)

public:
    /**
      * The algorithm approximates the BetweennessDir of all vertices so that the scores are
      * within an additive error @a epsilon with probability at least (1- @a delta).
      * The values are normalized by default.
      *
      * @param	G			the graph
      * @param  storePredecessors   keep track of the lists of predecessors?
      * @param	epsilon		maximum additive error
      * @param	delta		probability that the values are within the error guarantee
     */
    DynApproxBetweennessDir(const Graph& G, double epsilon=0.01, double delta=0.1, bool storePredecessors = true);

    /**
     * Runs the static approximated BetweennessDir centrality algorithm on the initial graph.
     */
    void run() override;

    /**
    * Updates the BetweennessDir centralities after a batch of edge insertions on the graph.
    *
    * @param batch The batch of edge insertions.
    */
    void update(const std::vector<GraphEvent>& batch);

    /**
    * Get number of path samples used for last calculation
    */
    count getNumberOfSamples();

    count computeVDdirected();

private:
    count DFS(node u, count j, Graph & sccDAG, std::vector<bool>& marked, std::vector<node>& sorted);
    void sampleNewPaths(count start, count end);
    bool storePreds = false; //TODO add vector with predecessors
    double epsilon; //!< maximum error
    double delta;
    const double c = 0.5; // universal positive constant - see reference in paper
    count r;
    count r2;
<<<<<<< HEAD
    std::vector<BFSvisit> sssp;
    std::vector<BFSvisit> sssp2;
=======
    std::vector<DynSSSP> sssp;
    std::vector<DynSSSP> sssp2;
>>>>>>> 7ca6f73b0... Fix includes and first steps at compiling files (not compiling yet)
    std::vector<std::vector<count>> npaths;

    std::vector<count> vis;
    std::vector<node> u;
    std::vector<node> v;
    std::vector<count> maxDist;
    std::vector<count> maxDist2;
    std::vector <std::vector<node>> sampledPaths;

    const count infDist = std::numeric_limits<count>::max();
};

} /* namespace NetworKit */

#endif // NETWORKIT_CENTRALITY_DYN_APPROX_BETWEENNESS_DIR_HPP_
