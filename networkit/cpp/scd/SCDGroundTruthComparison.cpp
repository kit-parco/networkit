#include <string>
#include <networkit/auxiliary/SignalHandling.hpp>
#include <networkit/scd/SCDGroundTruthComparison.hpp>

namespace NetworKit {
SCDGroundTruthComparison::SCDGroundTruthComparison(const Graph &g, const Cover &groundTruth,
                                                   const std::map<node, std::set<node>> &found,
                                                   bool ignoreSeeds)
    : g(&g), groundTruth(&groundTruth), found(&found), ignoreSeeds(ignoreSeeds) {}

void SCDGroundTruthComparison::run() {
    Aux::SignalHandler handler;

    hasRun = false;

    jaccardScores.clear();
    f1Scores.clear();
    precisionScores.clear();
    recallScores.clear();

    std::vector<std::vector<node>> groundTruthSets(groundTruth->upperBound());
    std::vector<count> foundSizes(g->upperNodeIdBound(), 0),
        truthSizes(groundTruth->upperBound(), 0);

    g->forNodes([&](node u) {
        for (index s : (*groundTruth)[u]) {
            groundTruthSets[s].emplace_back(u);
            ++truthSizes[s];
        }
    });

    for (const auto &it : *found) {
        for (node u : it.second) {
            if (g->hasNode(u)) {
                ++foundSizes[it.first];
            }
        }
    }

    handler.assureRunning();

    std::map<index, count> overlap;

    for (const auto &foundIt : *found) {
        handler.assureRunning();

        node seed = foundIt.first;

        if (!ignoreSeeds && !g->hasNode(seed)) {
            throw std::runtime_error("Error, the graph does not contain the seed node "
                                     + std::to_string(seed));
        }

        overlap.clear();

        const auto &foundNodes = foundIt.second;

        std::set<node> allowedSubsets;

        if (!ignoreSeeds) {
            allowedSubsets = groundTruth->subsetsOf(seed);
        }

        for (node u : foundNodes) {
            if (g->hasNode(u)) {
                for (index s : (*groundTruth)[u]) {
                    if (ignoreSeeds || allowedSubsets.count(s) > 0) {
                        ++overlap[s];
                    }
                }
            }
        }

        handler.assureRunning();

        double bestJaccard = 0;
        double bestF1 = 0;
        double bestPrecision = 0;
        double bestRecall = 0;

        for (auto o : overlap) {
            double currentj =
                static_cast<double>(o.second)
                / static_cast<double>(foundSizes[seed] + truthSizes[o.first] - overlap[o.first]);

            double recall = o.second * 1.0 / truthSizes[o.first];
            double precision = o.second * 1.0 / foundSizes[seed];

            double currentf1 = 0;
            if (precision > 0 && recall > 0) {
                currentf1 = 2 * (precision * recall) / (precision + recall);
            }

            if (currentj > bestJaccard) {
                bestJaccard = currentj;
            }

            if (currentf1 > bestF1) {
                bestF1 = currentf1;
            }

            if (recall > bestRecall) {
                bestRecall = recall;
            }

            if (precision > bestPrecision) {
                bestPrecision = precision;
            }
        }

        jaccardScores[seed] = bestJaccard;
        f1Scores[seed] = bestF1;
        precisionScores[seed] = bestPrecision;
        recallScores[seed] = bestRecall;
    }

    handler.assureRunning();

    auto avgScore = [](const std::map<index, double> &scores) -> double {
        double scoreSum = 0;
        for (auto a : scores) {
            scoreSum += a.second;
        }

        return scoreSum / scores.size();
    };

    averageJaccard = avgScore(jaccardScores);
    averageF1 = avgScore(f1Scores);
    averageRecall = avgScore(recallScores);
    averagePrecision = avgScore(precisionScores);

    handler.assureRunning();
    hasRun = true;
}

const std::map<index, double> &SCDGroundTruthComparison::getIndividualJaccard() const {
    assureFinished();
    return jaccardScores;
}

const std::map<index, double> &SCDGroundTruthComparison::getIndividualF1() const {
    assureFinished();
    return f1Scores;
}

const std::map<index, double> &SCDGroundTruthComparison::getIndividualRecall() const {
    assureFinished();
    return recallScores;
}

const std::map<index, double> &SCDGroundTruthComparison::getIndividualPrecision() const {
    assureFinished();
    return precisionScores;
}

double SCDGroundTruthComparison::getAverageJaccard() const {
    assureFinished();
    return averageJaccard;
}

double SCDGroundTruthComparison::getAverageF1() const {
    assureFinished();
    return averageF1;
}

double SCDGroundTruthComparison::getAverageRecall() const {
    assureFinished();
    return averageRecall;
}

double SCDGroundTruthComparison::getAveragePrecision() const {
    assureFinished();
    return averagePrecision;
}
} // namespace NetworKit
