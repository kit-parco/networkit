/*
 * EnsemblePreprocessing.cpp
 *
 *  Created on: 26.02.2013
 *      Author: Christian Staudt (christian.staudt@kit.edu)
 */

#include "EPP.h"

#include "../coarsening/ClusterContractor.h"
#include "../coarsening/ClusteringProjector.h"
#include "../community/JaccardMeasure.h"
#include "../auxiliary/Log.h"
#include "PLM.h"
#include "PLP.h"
#include "CNM.h"

namespace NetworKit {

EPP::EPP(const Graph& G) : CommunityDetectionAlgorithm(G) {
	this->finalClusterer = NULL;
	this->overlap = NULL;
}


void EPP::addBaseClusterer(std::unique_ptr<CommunityDetectionAlgorithm>& base) {
	this->baseClusterers.push_back(std::move(base));
}

void EPP::setFinalClusterer(std::unique_ptr<CommunityDetectionAlgorithm>& final) {
	this->finalClusterer = std::move(final);
}

void EPP::setOverlapper(std::unique_ptr<Overlapper>& overlap) {
	this->overlap = std::move(overlap);
}

void EPP::runImpl() {
	INFO("STARTING EnsemblePreprocessing on G=" , G.toString());

	// fixed sub-algorithms
	ClusterContractor contracter;
	ClusteringProjector projector;

	// data
	baseClusterings.clear();
	baseClusterings.resize(baseClusterers.size(), Partition(G.upperNodeIdBound())); // collection of base clusterings - fill with empty clustering

	// run base clusterers in parallel
	#pragma omp parallel for
	for (index b = 0; b < baseClusterers.size(); b += 1) {
		// FIXME: initialization of base clusterer?
		baseClusterers.at(b)->run();
		if (baseClusterers.at(b)->hasFinished())
			baseClusterings.at(b) = baseClusterers.at(b)->getPartition();
	}

	// ANALYSIS
	assureRunning();
	if (CALC_DISSIMILARITY) {
		JaccardMeasure dm;
		double dissimilaritySum = 0.0;
		for (index b = 0; b < baseClusterings.size(); b += 1) {
			for (index c = b + 1; c < baseClusterings.size(); c += 1) {
				double d = dm.getDissimilarity(G, baseClusterings.at(b), baseClusterings.at(c));
				dissimilaritySum += d;
			}
		}
		double avgDissimilarity = dissimilaritySum / (baseClusterings.size() * (baseClusterings.size() - 1) / 2.0);
		std::cout << "[INFO] avg. base clustering dissimilarity: " << avgDissimilarity << std::endl;
	}
	//

	// create core clustering
	core = this->overlap->run(G, baseClusterings);
	// contract graph according to core clustering
	std::pair<Graph, std::vector<node> > contraction = contracter.run(G, core);
	Graph Gcore = contraction.first;
	std::vector<node> fineToCoarse = contraction.second;
	// send contracted graph to final clusterer
	// TODO: maybe put this in a private helper function as this could be distracting...
	if (auto tmp = dynamic_cast<PLM*>(this->finalClusterer.get())) {
		DEBUG("final clusterer is PLM");
		this->finalClusterer.reset(new PLM(Gcore, tmp));
	} else if (auto tmp = dynamic_cast<PLP*>(this->finalClusterer.get())) {
		DEBUG("final clusterer is PLP");
		this->finalClusterer.reset(new PLP(Gcore, *tmp));
	} else if (auto tmp = dynamic_cast<CNM*>(this->finalClusterer.get())) {
		DEBUG("final clusterer is CNM");
		this->finalClusterer.reset(new CNM(Gcore));
	}
	assureRunning();
	this->finalClusterer->run();
	assureRunning();
		Partition finalCoarse = this->finalClusterer->getPartition();

		// project clustering of contracted graph back to original graph
		Partition final = projector.projectBack(Gcore, G, fineToCoarse, finalCoarse);
		// return clustering
		result = std::move(final);
}

std::string EPP::toString() const {
	std::stringstream strm;
	strm << "EnsemblePreprocessing(" << "base=" << this->baseClusterers.front()->toString() << ",ensemble=" << this->baseClusterers.size() << ",final=" << this->finalClusterer->toString() << ")";
	return strm.str();
}


Partition EPP::getCorePartition() const {
	return core;
}

std::vector<Partition> EPP::getBasePartitions() const {
	return baseClusterings;
}

} /* namespace NetworKit */
