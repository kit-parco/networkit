/*
 * EffectiveDiameter.cpp
 *
 *  Created on: 16.06.2014
 *      Author: Marc Nemes
 */

#include "EffectiveDiameter.h"
#include "ConnectedComponents.h"

#include <math.h>
#include <iterator>
#include <stdlib.h>
#include <omp.h>

namespace NetworKit {

double EffectiveDiameter::effectiveDiameter(const Graph& G) {
	return EffectiveDiameter::effectiveDiameter(G, 0.9, 64, 7);
}

double EffectiveDiameter::effectiveDiameter(const Graph& G, const double ratio) {
	return EffectiveDiameter::effectiveDiameter(G, ratio, 64, 7);
}

/*
this is a variaton of the ANF algorithm presented in the paper "A Fast and Scalable Tool for Data Mining
in Massive Graphs" by Palmer, Gibbons and Faloutsos which can be found here: http://www.cs.cmu.edu/~christos/PUBLICATIONS/kdd02-anf.pdf
*/
double EffectiveDiameter::effectiveDiameter(const Graph& G, const double ratio, const count k, const count r) {
	// the length of the bitmask where the number of connected nodes is saved
	count lengthOfBitmask = (count) ceil(log2(G.numberOfNodes()));
	// saves all k bitmasks for every node of the current iteration
	std::vector<std::vector<unsigned int> > mCurr;
	// saves all k bitmasks for every node of the previous iteration
	std::vector<std::vector<unsigned int> > mPrev;
	// the list of nodes that are already connected to all other nodes
	std::vector<node> finishedNodes;
	// the maximum possible bitmask based on the random initialization of all k bitmasks
	std::vector<count> highestCount;
	// the amount of nodes that need to be connected to all others nodes
	count threshold = (count) (ceil(ratio * G.numberOfNodes()));
	// the current distance of the neighborhoods
	count h = 1;
	// the amount of nodes that are connected to all other nodes (|finishedNodes|)
	count numberOfFinishedNodes = 0;
	// sums over the number of edges needed to reach 90% of all other nodes
	double effectiveDiameter = 0;
	// the estimated number of connected nodes
	double estimatedConnectedNodes;

	double random;
	srand (time(NULL));

	// initialize all vectors
	for (count j = 0; j < k; j++) {
		highestCount.push_back(j);
		highestCount[j] = 0;
	}
	G.forNodes([&](node v) {
		finishedNodes.push_back(v);
		finishedNodes[v] = 0;
		std::vector<unsigned int> tmp;
		for (count j = 0; j < k; j++) {
			tmp.push_back(0);
		}
		mCurr.push_back(tmp);
		mPrev.push_back(tmp);

		// set one bit in each bitmask with probability P(bit i=1) = 0.5^(i+1), i=0,..
		for (count j = 0; j < k; j++) {
			random = (rand()/(double)(RAND_MAX));
			for (count i = 0; i < lengthOfBitmask+r; i++) {
				if (random > pow(0.5,i+1)) {
					mPrev[v][j] |= 1 << i;
					break;
				}
			}
			// add the current bit to the maximum-bitmask
			highestCount[j] = highestCount[j] | mPrev[v][j];
		}
	});

	// as long as we need to connect more nodes
	while (numberOfFinishedNodes < G.numberOfNodes()) {
		G.forNodes([&](node v) {
			// if the current node is not yet connected to all other nodes
			if (finishedNodes[v] == 0) {
				#pragma omp parallel for
				// for each parallel approximation
				for (count j = 0; j < k; j++) {
					// the node is still connected to all previous neighbors
					mCurr[v][j] = mPrev[v][j];
					// and to all previous neighbors of all its neighbors
					G.forNeighborsOf(v, [&](node u) {
						mCurr[v][j] = mCurr[v][j] | mPrev[u][j];
					});
				}
				// the least bit number in the bitmask of the current node/distance that has not been set
				double b = 0;

				for (count j = 0; j < k; j++) {
					for (count i = 0; i < sizeof(i)*8; i++) {
						if (((mCurr[v][j] >> i) & 1) == 0) {
							b += i;
							break;
						}
					}
				}
				// calculate the average least bit number that has not been set over all parallel approximations
				b = b / k;
				// calculate the estimated number of neighbors
				estimatedConnectedNodes = (pow(2,b) / 0.77351);

				// check whether all k bitmask for this node have reached their highest possible value
				bool nodeFinished = true;
				for (count j = 0; j < k; j++) {
					if (mCurr[v][j] != highestCount[j]) {
						nodeFinished = false;
						break;
					}
				}
				// if the node wont change or is connected to enough nodes it must no longer be considered
				if (estimatedConnectedNodes >= threshold || nodeFinished) {
					finishedNodes[v] = 1;
					numberOfFinishedNodes++;
					effectiveDiameter += h;
				}
			}
		});
		mPrev = mCurr;
		h++;
	}
	return effectiveDiameter/G.numberOfNodes();
}

double EffectiveDiameter::effectiveDiameterExact(const Graph& G) {
	return EffectiveDiameter::effectiveDiameterExact(G, 0.9);
}

/*
this is a variaton of the ANF algorithm presented in the paper "A Fast and Scalable Tool for Data Mining
in Massive Graphs" by Palmer, Gibbons and Faloutsos which can be found here: http://www.cs.cmu.edu/~christos/PUBLICATIONS/kdd02-anf.pdf
*/
double EffectiveDiameter::effectiveDiameterExact(const Graph& G, const double ratio) {
	// list of nodes that are already connected to enough other nodes
	std::vector<node> finishedNodes;
	// saves the reachable nodes of the current iteration
	std::vector<std::vector<bool> > mCurr;
	// saves the reachable nodes of the previous iteration
	std::vector<std::vector<bool> > mPrev;
	// the amount of nodes that are connected to enough other nodes (|finishedNodes|)
	count numberOfFinishedNodes = 0;
	// sums over the number of edges needed to reach 90% of all other nodes
	double effectiveDiameter = 0;
	// the current distance of the neighborhoods
	count h = 1;
	// number of nodes that need to be connected with all other nodes
	count threshold = (uint64_t) (ceil(ratio * G.numberOfNodes()) + 0.5);

	// initialize all nodes
	G.forNodes([&](node v){
		std::vector<bool> connectedNodes;
		// initialize n entries with value 0
		connectedNodes.assign(G.numberOfNodes(),0);
		// the node is always connected to itself
		connectedNodes[v] = 1;
		finishedNodes.push_back(v);
		finishedNodes[v] = 0;
		mCurr.push_back(connectedNodes);
		mPrev.push_back(connectedNodes);
	});

	// as long as we need to connect more nodes
	while (numberOfFinishedNodes < G.numberOfNodes()) {
			G.forNodes([&](node v) {
				// if the current node is not yet connected to all other nodes
				if (finishedNodes[v] == 0) {
					mCurr[v] = mPrev[v];
					G.forNeighborsOf(v, [&](node u) {
						for (count i = 0; i < G.numberOfNodes(); i++) {
							// add the current neighbor of u to the neighborhood of v
							mCurr[v][i] = mCurr[v][i] || mPrev[u][i];
						}
					});

					// compute the number of connected nodes
					count numConnectedNodes = 0;
					for (count i = 0; i < G.numberOfNodes(); i++) {
						if (mCurr[v][i] == 1) {
							numConnectedNodes++;
						}
					}

					// when the number of connected nodes surpasses the threshold the node must no longer be considered
					if (numConnectedNodes >= threshold) {
						finishedNodes[v] = 1;
						numberOfFinishedNodes++;
						effectiveDiameter += h;
					}
				}
			});
			mPrev = mCurr;
			h++;
		}
	// return the found effective diameter
	return effectiveDiameter/G.numberOfNodes();
}

/*
this is a variaton of the ANF algorithm presented in the paper "A Fast and Scalable Tool for Data Mining
in Massive Graphs" by Palmer, Gibbons and Faloutsos which can be found here: http://www.cs.cmu.edu/~christos/PUBLICATIONS/kdd02-anf.pdf
*/
std::map<count, double> EffectiveDiameter::hopPlot(const Graph& G, const count maxDistance, const count k, const count r) {
	//the returned hop-plot
	std::map<count, double> hopPlot;
	// the length of the bitmask where the number of connected nodes is saved
	count lengthOfBitmask = (count) ceil(log2(G.numberOfNodes()));
	// saves all k bitmasks for every node of the current iteration
	std::vector<std::vector<unsigned int> > mCurr;
	// saves all k bitmasks for every node of the previous iteration
	std::vector<std::vector<unsigned int> > mPrev;
	// the list of nodes that are already connected to all other nodes
	std::vector<node> finishedNodes;
	// the maximum possible bitmask based on the random initialization of all k bitmasks
	std::vector<count> highestCount;
	// the current distance of the neighborhoods
	count h = 1;
	// the amount of nodes that are connected to all other nodes (|finishedNodes|)
	count numberOfFinishedNodes = 0;
	// the estimated number of connected nodes
	double estimatedConnectedNodes;
	// the sum over all estimated connected nodes
	double totalConnectedNodes;

	double random;
	srand (time(NULL));

	// initialize all vectors
	for (count j = 0; j < k; j++) {
		highestCount.push_back(j);
		highestCount[j] = 0;
	}
	G.forNodes([&](node v) {
		finishedNodes.push_back(v);
		finishedNodes[v] = 0;
		std::vector<unsigned int> tmp;
		for (count j = 0; j < k; j++) {
			tmp.push_back(0);
		}
		mCurr.push_back(tmp);
		mPrev.push_back(tmp);

		// set one bit in each bitmask with probability P(bit i=1) = 0.5^(i+1), i=0,..
		for (count j = 0; j < k; j++) {
			random = (rand()/(double)(RAND_MAX));
			for (count i = 0; i < lengthOfBitmask+r; i++) {
				if (random > pow(0.5,i+1)) {
					mPrev[v][j] |= 1 << i;
					break;
				}
			}
			// add the current bit to the maximum-bitmask
			highestCount[j] = highestCount[j] | mPrev[v][j];
		}
	});
	// at zero distance, all nodes can only reach themselves
	hopPlot[0] = 1/G.numberOfNodes();

	// as long as we need to connect more nodes
	while (numberOfFinishedNodes < G.numberOfNodes() && (maxDistance <= 0 || h < maxDistance)) {
		totalConnectedNodes = 0;
		G.forNodes([&](node v) {
			// if the current node is not yet connected to all other nodes
			if (finishedNodes[v] == 0) {
				// for each parallel approximation
				for (count j = 0; j < k; j++) {
					// the node is still connected to all previous neighbors
					mCurr[v][j] = mPrev[v][j];
					// and to all previous neighbors of all its neighbors
					G.forNeighborsOf(v, [&](node u) {
						mCurr[v][j] = mCurr[v][j] | mPrev[u][j];
					});
				}
				// the least bit number in the bitmask of the current node/distance that has not been set
				double b = 0;

				for (count j = 0; j < k; j++) {
					for (count i = 0; i < sizeof(i)*8; i++) {
						if (((mCurr[v][j] >> i) & 1) == 0) {
							b += i;
							break;
						}
					}
				}
				// calculate the average least bit number that has not been set over all parallel approximations
				b = b / k;
				// calculate the estimated number of neighbors
				estimatedConnectedNodes = (pow(2,b) / 0.77351);

				// enforce monotonicity
				if (estimatedConnectedNodes > G.numberOfNodes()) {
					estimatedConnectedNodes = G.numberOfNodes();
				}

				// add value of the node to all nodes so we can calculate the average
				totalConnectedNodes += estimatedConnectedNodes;

				// check whether all k bitmask for this node have reached the highest possible value
				bool nodeFinished = true;
				for (count j = 0; j < k; j++) {
					if (mCurr[v][j] != highestCount[j]) {
						nodeFinished = false;
						break;
					}
				}

				// if the node wont change or is connected to enough nodes it must no longer be considered
				if (estimatedConnectedNodes >= G.numberOfNodes() || nodeFinished) {
					finishedNodes[v] = 1;
					numberOfFinishedNodes++;
				}
			} else {
				totalConnectedNodes += G.numberOfNodes();
			}
		});
		hopPlot[h] = totalConnectedNodes/(G.numberOfNodes()*G.numberOfNodes());
		mPrev = mCurr;
		h++;
	}
	return hopPlot;
}

}
