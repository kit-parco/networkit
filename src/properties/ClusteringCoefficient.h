/*
 * ClusteringCoefficient.h
 *
 *  Created on: 08.04.2013
 *      Author: Lukas Barth, David Weiß
 */

#ifndef CLUSTERINGCOEFFICIENT_H_
#define CLUSTERINGCOEFFICIENT_H_

#include "../graph/Graph.h"

namespace NetworKit {


class ClusteringCoefficient {

public:

	ClusteringCoefficient();

	virtual ~ClusteringCoefficient();

	std::vector<double> local(Graph &G) const;
	
	/**
	 * This calculates the average local clustering coefficient
	 * $$c(G) := \frac{1}{n} \sum_{u \in V} c(u)$$
	 *
	 * where $c(u) := \frac{2 \cdot |E(N(u))| }{\deg(u) \cdot ( \deg(u) - 1)}$
	 */
	virtual double avgLocal(Graph& G) const;
  	virtual double approxAvgLocal(Graph& G, const count tries) const;

	/**
	 * This calculates the global clustering coefficient
	 */
  	virtual double global(Graph& G) const;
  	virtual double approxGlobal(Graph& G, const count tries) const;
  
};

} /* namespace NetworKit */
#endif /* CLUSTERINGCOEFFICIENT_H_ */
