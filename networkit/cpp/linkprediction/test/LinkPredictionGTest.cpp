// no-networkit-format
/*
 * LinkPredictionGTest.cpp
 *
 *  Created on: 07.12.2014
 *      Author: Kolja Esders (kolja.esders@student.kit.edu)
 */

#include <gtest/gtest.h>

#include <chrono>

#include <networkit/graph/Graph.hpp>

#include <networkit/io/METISGraphReader.hpp>
#include <networkit/linkprediction/KatzIndex.hpp>
#include <networkit/linkprediction/CommonNeighborsIndex.hpp>
#include <networkit/linkprediction/JaccardIndex.hpp>
#include <networkit/linkprediction/ROCMetric.hpp>
#include <networkit/linkprediction/PrecisionRecallMetric.hpp>
#include <networkit/linkprediction/RandomLinkSampler.hpp>
#include <networkit/linkprediction/MissingLinksFinder.hpp>
#include <networkit/linkprediction/UDegreeIndex.hpp>
#include <networkit/linkprediction/VDegreeIndex.hpp>
#include <networkit/linkprediction/LinkThresholder.hpp>
#include <networkit/linkprediction/TotalNeighborsIndex.hpp>
#include <networkit/linkprediction/NeighborsMeasureIndex.hpp>
#include <networkit/linkprediction/SameCommunityIndex.hpp>
#include <networkit/linkprediction/PredictionsSorter.hpp>
#include <networkit/linkprediction/LinkPredictor.hpp>

namespace NetworKit {

class LinkPredictionGTest : public testing::Test {
public:
    void SetUp();

protected:
    Graph G{7};
    Graph trainingGraph;

    std::vector<std::pair<node, node>> missingLinks;
    std::vector<LinkPredictor::prediction> predictions;
};

void LinkPredictionGTest::SetUp() {
  G.addEdge(0, 1);
  G.addEdge(0, 3);
  G.addEdge(1, 2);
  G.addEdge(1, 4);
  G.addEdge(2, 3);
  G.addEdge(2, 4);
  G.addEdge(2, 5);
  G.addEdge(3, 4);
  G.addEdge(3, 5);
  G.addEdge(4, 5);
  trainingGraph = G;
  trainingGraph.removeEdge(0, 1);
  trainingGraph.removeEdge(2, 4);
  trainingGraph.removeEdge(3, 5);
  missingLinks = MissingLinksFinder(trainingGraph).findAtDistance(2);
  CommonNeighborsIndex cn(trainingGraph);
  predictions = cn.runOn(missingLinks);
  PredictionsSorter::sortByScore(predictions);
}

TEST_F(LinkPredictionGTest, testCommonNeighborsIndexRunOn) {
  EXPECT_EQ(6, predictions.size());
  EXPECT_EQ(2, predictions[0].first.first); EXPECT_EQ(4, predictions[0].first.second); EXPECT_EQ(3, predictions[0].second);
  EXPECT_EQ(1, predictions[1].first.first); EXPECT_EQ(3, predictions[1].first.second); EXPECT_EQ(2, predictions[1].second);
  EXPECT_EQ(1, predictions[2].first.first); EXPECT_EQ(5, predictions[2].first.second); EXPECT_EQ(2, predictions[2].second);
  EXPECT_EQ(3, predictions[3].first.first); EXPECT_EQ(5, predictions[3].first.second); EXPECT_EQ(2, predictions[3].second);
  EXPECT_EQ(0, predictions[4].first.first); EXPECT_EQ(2, predictions[4].first.second); EXPECT_EQ(1, predictions[4].second);
  EXPECT_EQ(0, predictions[5].first.first); EXPECT_EQ(4, predictions[5].first.second); EXPECT_EQ(1, predictions[5].second);
}

TEST_F(LinkPredictionGTest, testMissingLinksFinderDistanceTwo) {
  EXPECT_EQ(6, missingLinks.size());
  EXPECT_EQ(0, missingLinks[0].first); EXPECT_EQ(2, missingLinks[0].second);
  EXPECT_EQ(0, missingLinks[1].first); EXPECT_EQ(4, missingLinks[1].second);
  EXPECT_EQ(1, missingLinks[2].first); EXPECT_EQ(3, missingLinks[2].second);
  EXPECT_EQ(1, missingLinks[3].first); EXPECT_EQ(5, missingLinks[3].second);
  EXPECT_EQ(2, missingLinks[4].first); EXPECT_EQ(4, missingLinks[4].second);
  EXPECT_EQ(3, missingLinks[5].first); EXPECT_EQ(5, missingLinks[5].second);
}

TEST_F(LinkPredictionGTest, testMissingLinksFinderDistanceThree) {
  std::vector<std::pair<node, node>> hopThreeMissingLinks = MissingLinksFinder(trainingGraph).findAtDistance(3);
  EXPECT_EQ(2, hopThreeMissingLinks.size());
  EXPECT_EQ(0, hopThreeMissingLinks[0].first); EXPECT_EQ(1, hopThreeMissingLinks[0].second);
  EXPECT_EQ(0, hopThreeMissingLinks[1].first); EXPECT_EQ(5, hopThreeMissingLinks[1].second);
}

TEST_F(LinkPredictionGTest, testLinkThresholderByScore) {
  std::vector<std::pair<node, node>> selectedLinks;
  selectedLinks = LinkThresholder::byScore(predictions, 2);
  EXPECT_EQ(4, selectedLinks.size());
  EXPECT_EQ(1, selectedLinks[0].first); EXPECT_EQ(3, selectedLinks[0].second);
  EXPECT_EQ(1, selectedLinks[1].first); EXPECT_EQ(5, selectedLinks[1].second);
  EXPECT_EQ(2, selectedLinks[2].first); EXPECT_EQ(4, selectedLinks[2].second);
  EXPECT_EQ(3, selectedLinks[3].first); EXPECT_EQ(5, selectedLinks[3].second);
}

TEST_F(LinkPredictionGTest, testLinkThresholderByCount) {
  std::vector<std::pair<node, node>> selectedLinks;
  selectedLinks = LinkThresholder::byCount(predictions, 5);
  EXPECT_EQ(5, selectedLinks.size());
  EXPECT_EQ(0, selectedLinks[0].first); EXPECT_EQ(2, selectedLinks[0].second);
  EXPECT_EQ(1, selectedLinks[1].first); EXPECT_EQ(3, selectedLinks[1].second);
  EXPECT_EQ(1, selectedLinks[2].first); EXPECT_EQ(5, selectedLinks[2].second);
  EXPECT_EQ(2, selectedLinks[3].first); EXPECT_EQ(4, selectedLinks[3].second);
  EXPECT_EQ(3, selectedLinks[4].first); EXPECT_EQ(5, selectedLinks[4].second);
}

TEST_F(LinkPredictionGTest, testLinkThresholderByPercentage) {
  std::vector<std::pair<node, node>> selectedLinks;
  selectedLinks = LinkThresholder::byPercentage(predictions, 0.5);
  EXPECT_EQ(3, selectedLinks.size());
  EXPECT_EQ(1, selectedLinks[0].first); EXPECT_EQ(3, selectedLinks[0].second);
  EXPECT_EQ(1, selectedLinks[1].first); EXPECT_EQ(5, selectedLinks[1].second);
  EXPECT_EQ(2, selectedLinks[2].first); EXPECT_EQ(4, selectedLinks[2].second);
}

TEST_F(LinkPredictionGTest, testTrainingGraphGenerator) {
  Graph trainingGraph = RandomLinkSampler::byPercentage(G, 0.7);
  EXPECT_EQ(7, trainingGraph.numberOfEdges());
}

TEST_F(LinkPredictionGTest, testUDegreeIndexRun) {
  UDegreeIndex uDegreeIndex(G);
  EXPECT_EQ(2, uDegreeIndex.run(0, 3));
  EXPECT_EQ(3, uDegreeIndex.run(1, 0));
  EXPECT_EQ(4, uDegreeIndex.run(2, 3));
  EXPECT_EQ(4, uDegreeIndex.run(3, 4));
  EXPECT_EQ(4, uDegreeIndex.run(4, 5));
  EXPECT_EQ(3, uDegreeIndex.run(5, 4));
}

TEST_F(LinkPredictionGTest, testVDegreeIndexRun) {
  VDegreeIndex vDegreeIndex(G);
  EXPECT_EQ(2, vDegreeIndex.run(3, 0));
  EXPECT_EQ(3, vDegreeIndex.run(0, 1));
  EXPECT_EQ(4, vDegreeIndex.run(3, 2));
  EXPECT_EQ(4, vDegreeIndex.run(4, 3));
  EXPECT_EQ(4, vDegreeIndex.run(5, 4));
  EXPECT_EQ(3, vDegreeIndex.run(4, 5));
}

TEST_F(LinkPredictionGTest, testTotalNeighborsIndexRun) {
  TotalNeighborsIndex totalNeighborsIndex(trainingGraph);
  EXPECT_EQ(3, totalNeighborsIndex.run(0, 2));
  EXPECT_EQ(3, totalNeighborsIndex.run(0, 4));
  EXPECT_EQ(3, totalNeighborsIndex.run(1, 3));
  EXPECT_EQ(2, totalNeighborsIndex.run(1, 5));
  EXPECT_EQ(3, totalNeighborsIndex.run(2, 4));
  EXPECT_EQ(3, totalNeighborsIndex.run(3, 5));
}

TEST_F(LinkPredictionGTest, testNeighborsMeasureIndexRun) {
  NeighborsMeasureIndex neighborsMeasureIndex(trainingGraph);
  EXPECT_EQ(1, neighborsMeasureIndex.run(0, 2));
  EXPECT_EQ(1, neighborsMeasureIndex.run(0, 4));
  EXPECT_EQ(2, neighborsMeasureIndex.run(1, 3));
  EXPECT_EQ(2, neighborsMeasureIndex.run(1, 5));
  EXPECT_EQ(3, neighborsMeasureIndex.run(2, 4));
  EXPECT_EQ(2, neighborsMeasureIndex.run(3, 5));
}

TEST_F(LinkPredictionGTest, testROCMetric) {
  ROCMetric roc(G);
  std::pair<std::vector<double>, std::vector<double>> curve = roc.getCurve(predictions);
  double auc = roc.getAreaUnderCurve();
  EXPECT_EQ(auc, 0.8125);
  EXPECT_EQ(0, curve.first[0]); EXPECT_EQ(0.5, curve.second[0]);
  EXPECT_EQ(0.25, curve.first[1]); EXPECT_EQ(0.5, curve.second[1]);
  EXPECT_EQ(0.5, curve.first[2]); EXPECT_EQ(1, curve.second[2]);
  EXPECT_EQ(0.75, curve.first[3]); EXPECT_EQ(1, curve.second[3]);
  EXPECT_EQ(1, curve.first[4]); EXPECT_EQ(1, curve.second[4]);
}

TEST_F(LinkPredictionGTest, testPRMetric) {
  PrecisionRecallMetric pr(G);
  std::pair<std::vector<double>, std::vector<double>> curve = pr.getCurve(predictions);
  double auc = pr.getAreaUnderCurve();
  EXPECT_EQ(auc, 0.5);
  EXPECT_EQ(0, curve.first[0]); EXPECT_EQ(1, curve.second[0]);
  EXPECT_EQ(0.5, curve.first[1]); EXPECT_EQ(1.0 / 3, curve.second[1]);
  EXPECT_EQ(1, curve.first[2]); EXPECT_EQ(1.0 / 3, curve.second[2]);
}

TEST_F(LinkPredictionGTest, testKatzRunOnOrdering) {
  METISGraphReader graphReader;
  Graph newG = graphReader.read("input/jazz.graph");
  KatzIndex katz(newG);
  Graph trainingGraph = RandomLinkSampler::byPercentage(newG, 0.7);
  std::vector<std::pair<node, node>> nodePairs = MissingLinksFinder(trainingGraph).findAtDistance(2);
  std::vector<std::pair<std::pair<node, node>, double>> preds = katz.runOn(missingLinks);
  for (index i = 0; i < preds.size() - 1; ++i) {
    EXPECT_TRUE(preds[i].first.first < preds[i+1].first.first || preds[i].first.second < preds[i+1].first.second);
  }
}

} // namespace NetworKit
