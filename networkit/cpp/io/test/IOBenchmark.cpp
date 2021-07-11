// no-networkit-format
/*
 * IOBenchmark.cpp
 *
 *  Created on: 01.02.2013
 *      Author: Christian Staudt (christian.staudt@kit.edu)
 */

#include <gtest/gtest.h>
#include <vector>
#include <string>
#include <utility>
#include <algorithm>
#include <fstream>

#include <networkit/auxiliary/Log.hpp>
#include <networkit/auxiliary/Timer.hpp>
#include <networkit/io/METISGraphReader.hpp>

#include <networkit/io/RasterReader.hpp>
#include <networkit/generators/quadtree/QuadtreePolarEuclid.hpp>
#include <networkit/generators/quadtree/QuadtreeCartesianEuclid.hpp>
#include <networkit/geometric/HyperbolicSpace.hpp>

namespace NetworKit {

class IOBenchmark: public testing::Test {
public:
    static void convertToHeatMap(std::vector<bool> &infected, std::vector<double> &xcoords,
                                 std::vector<double> &ycoords, const std::string &filename,
                                 double resolution = 1);
};

/**
 * This should actually not be in the benchmark but somewhere else. Can't figure out where yet.
 */
void IOBenchmark::convertToHeatMap(std::vector<bool> &infected, std::vector<double> &xcoords,
                                   std::vector<double> &ycoords, const std::string &filename,
                                   double resolution) {

    auto minmaxx = std::minmax_element (xcoords.begin(),xcoords.end());
    auto  minmaxy = std::minmax_element (ycoords.begin(),ycoords.end());

    //create data structures for heat map
    double xspread = *minmaxx.second - *minmaxx.first;
    double yspread = *minmaxy.second - *minmaxy.first;
    int xsteps = xspread / resolution + 1;
    int ysteps = yspread / resolution + 1;
    std::vector<std::vector<int> > infectedByRegion(xsteps);
    std::vector<std::vector<bool> > populated(xsteps);
    for (int i = 0; i < xsteps; i++) {
        infectedByRegion[i].resize(ysteps, 0);
        populated[i].resize(ysteps, false);
    }

    //fill heat map
    for (index i = 0; i < infected.size(); i++) {
        //get bins
        index xBin = (xcoords[i]-*minmaxx.first)/resolution;
        index yBin = (ycoords[i]-*minmaxy.first)/resolution;

        //check bins
        assert(xBin < infectedByRegion.size());
        assert(yBin < infectedByRegion[0].size());
        populated[xBin][yBin] = true;

        if (infected[i]) {
            //write in bin
            infectedByRegion[xBin][yBin]++;
        }
    }

    //open output file
    std::ofstream file;
    file.open(filename.c_str());

    //write column labels
    file << "x" << "\t" << "y" << "\t" << "label" << std::endl;

    //write heat map
    for (index i = 0; i < infectedByRegion.size(); i++) {
        double x = *minmaxx.first + i*resolution;
        for (index j = 0; j < infectedByRegion[i].size(); j++) {
            double y = *minmaxy.first + j*resolution;
            std::string heat = populated[i][j] ? std::to_string(infectedByRegion[i][j]) : std::string("NaN");
            file << x << '\t' << y << '\t' << heat << std::endl;
        }
        file << std::endl;
    }
    file.close();
}

TEST_F(IOBenchmark, timeMETISGraphReader) {
    std::string path = "input/caidaRouterLevel.graph";
    Aux::Timer runtime;

    INFO("[BEGIN] reading graph: " , path);
    runtime.start();
    METISGraphReader reader;
    Graph G = reader.read(path);
    runtime.stop();
    INFO("[DONE] reading graph " , runtime.elapsedTag());

    EXPECT_TRUE(! G.isEmpty());

}


TEST_F(IOBenchmark, benchRasterReader) {
    double normalizationFactor = 0.05;
    RasterReader reader(normalizationFactor);
    std::vector<double> xcoords;
    std::vector<double> ycoords;
    Aux::Timer runtime;

    std::vector<std::string> countries = {"deu", "fra", "usa"};
    for (auto country: countries) {
        std::string path("input/" + country + "p00ag.asc");

        // read raster file
        INFO("[BEGIN] reading raster data set: ", path);
        runtime.start();
        std::tie(xcoords, ycoords) = reader.read(path);
        runtime.stop();
        INFO("[DONE] reading raster data set with ", xcoords.size(), " points:" , runtime.elapsedTag());
        EXPECT_EQ(xcoords.size(), ycoords.size());

        const count n = xcoords.size();
        count numRuns = 1;

        for (index run = 0; run < numRuns; run++) {
            //transform into polar coordinates
            runtime.start();
            std::vector<double> angles(n);
            std::vector<double> radii(n);
            std::vector<Point2DWithIndex<double> > positions(n);
            std::vector<index> content(n);
            double maxR = 0;
            for (index i = 0; i < n; i++) {
                HyperbolicSpace::cartesianToPolar(Point2DWithIndex<double>(xcoords[i], ycoords[i]), angles[i], radii[i]);
                maxR = std::max(maxR, radii[i]);
                positions[i] = Point2DWithIndex<double>(xcoords[i], ycoords[i]);
                content[i] = i;
            }
            runtime.stop();

            INFO("Converted coordinates", runtime.elapsedTag());
            //define query function
            auto edgeProb = [n](double distance) -> double {return (1/distance)*exp(5)/(double)n ;};
            //auto edgeProb = [beta, thresholdDistance](double distance) -> double {return 1 / (exp(beta*(distance-thresholdDistance)/2)+1);};

            // construct quadtree
            runtime.start();
            QuadtreePolarEuclid<index> tree(angles, radii, content);
            runtime.stop();
            INFO("Filled quadtree", runtime.elapsedTag());

            uint64_t numQueries = 1000;
            long treeTotalNeighbours = 0;

            // perform range queries
            runtime.start();
            for (uint64_t q = 0; q < numQueries; ++q) {
                std::vector<index> result;
                index comparison = Aux::Random::integer(xcoords.size());
                Point2DWithIndex<double> query(xcoords[comparison], ycoords[comparison]);
                tree.getElementsProbabilistically(query, edgeProb, result);
                treeTotalNeighbours += result.size();
            }
            runtime.stop();
            INFO("Completed ", numQueries, " quadtree queries with on average ", treeTotalNeighbours / numQueries, " neighbours", runtime.elapsedTag());

            long naiveTotalNeighbours = 0;
            runtime.start();
            for (uint64_t q = 0; q < numQueries; ++q) {
                std::vector<index> result;
                index comparison = Aux::Random::integer(xcoords.size());
                double x = xcoords[comparison];
                double y = ycoords[comparison];
                for (index i = 0; i < xcoords.size(); i++) {
                    double xdiff = xcoords[i] - x;
                    double ydiff = ycoords[i] - y;
                    double prob = edgeProb(pow(xdiff*xdiff+ydiff*ydiff , 0.5));
                    double random = Aux::Random::real();
                    if (random < prob) result.push_back(i);
                }
                naiveTotalNeighbours += result.size();
            }
            runtime.stop();
            INFO("Completed ", numQueries, " naive queries with on average ", naiveTotalNeighbours / numQueries, " neighbours", runtime.elapsedTag());
            EXPECT_NEAR(treeTotalNeighbours, naiveTotalNeighbours, treeTotalNeighbours / 10);
        }
    }
}

/**
 * This test simulates a disease progression with the probabilistic range query. We use the susceptible-infected-susceptible model.
 */
TEST_F(IOBenchmark, simulateDiseaseProgression) {
    const double recoveryProb = 0.8;
    double normalizationFactor = 0.05;
    RasterReader reader(normalizationFactor);
    std::vector<double> xcoords;
    std::vector<double> ycoords;
    Aux::Timer runtime;

    std::vector<std::string> countries = {"deu", "fra", "usa"};
    for (auto country: countries) {
        std::string path("input/" + country + "p00ag.asc");

        // read raster file
        INFO("[BEGIN] reading raster data set: ", path);
        runtime.start();
        std::tie(xcoords, ycoords) = reader.read(path);
        runtime.stop();
        INFO("[DONE] reading raster data set " , runtime.elapsedTag());
        EXPECT_EQ(xcoords.size(), ycoords.size());
        const count n = xcoords.size();

        // get minimum and maximum
        auto minmaxx = std::minmax_element (xcoords.begin(),xcoords.end());
        INFO("X coordinates range from ", *minmaxx.first, " to ", *minmaxx.second, ".");
        auto  minmaxy = std::minmax_element (ycoords.begin(),ycoords.end());
        INFO("Y coordinates range from ", *minmaxy.first, " to ", *minmaxy.second, ".");

        //set resolution for heat map output
        double resolution;
        if (country.compare("usa") == 0) {
            resolution = 20;
        } else {
            resolution = 3;
        }

        //set neighbor probability
        auto edgeProb = [](double distance) -> double {return (1/distance)*exp(-11);};

        //convert coordinates
        runtime.start();
        std::vector<double> angles(n);
        std::vector<double> radii(n);
        std::vector<index> content(n);
        double maxR = 0;

        for (index i = 0; i < n; i++) {
            HyperbolicSpace::cartesianToPolar(Point2DWithIndex<double>(xcoords[i], ycoords[i]), angles[i], radii[i]);
            maxR = std::max(maxR, radii[i]);
            content[i] = i;
        }
        runtime.stop();

        // construct quadtree
        runtime.start();
        QuadtreePolarEuclid<index> tree(angles, radii, content);
        runtime.stop();
        INFO("Filled quadtree", runtime.elapsedTag());

        //data structures
        std::vector<bool> wasEverInfected(n, false);
        std::vector<bool> infectedState(n, false);
        std::vector<index> infectedList;

        //patient zero
        runtime.start();
        index patientZero = Aux::Random::index(n);
        INFO("Person ", patientZero, " at (", xcoords[patientZero], ", ", ycoords[patientZero], ") selected as patient zero.");
        infectedList.push_back(patientZero);
        infectedState[patientZero] = true;
        wasEverInfected[patientZero] = true;
        index step = 0;

        //main loop
        while (!infectedList.empty()) {
            INFO("At step ", step, ", ", infectedList.size(), " people are infected.");

            //get new infections
            std::vector<index> newInfections;
            for (index patient : infectedList) {
                Point2DWithIndex<double> query(xcoords[patient], ycoords[patient]);
                tree.getElementsProbabilistically(query, edgeProb, newInfections);
            }

            //remove old infections from list of new infections
            std::sort(newInfections.begin(), newInfections.end());
            newInfections.erase(std::unique(newInfections.begin(), newInfections.end()), newInfections.end());
            count newInfectionCount = newInfections.size();
            newInfections.erase(std::remove_if(newInfections.begin(), newInfections.end(), [wasEverInfected](index i)->bool{return wasEverInfected[i];}), newInfections.end());
            INFO(newInfectionCount, " infections happened, of which ", newInfectionCount - newInfections.size(), " were already infected or immune.");

            //get coordinates of new infections
            if (newInfections.size() > 0) {
                double minX = xcoords[newInfections[0]];
                double maxX = xcoords[newInfections[0]];
                double minY = ycoords[newInfections[0]];
                double maxY = ycoords[newInfections[0]];
                for (index patient : newInfections) {
                    minX = min(xcoords[patient], minX);
                    maxX = max(xcoords[patient], maxX);
                    minY = min(ycoords[patient], minY);
                    maxY = max(ycoords[patient], maxY);
                }
                INFO("X coordinates of new infections range from ", minX, " to ", maxX, ".");
                INFO("Y coordinates of new infections range from ", minY, " to ", maxY, ".");
            }

            //old infections may recover or stay infectious
            for (index oldPatient : infectedList) {
                if (Aux::Random::real() < recoveryProb) {
                    //patient has recovered
                    infectedState[oldPatient] = false;
                } else {
                    //patient still sick. Since the previous patients had been removed from the list of new infections, there won't be any duplicates.
                    newInfections.push_back(oldPatient);
                }
            }

            infectedList = newInfections;

            for (index patient : infectedList) {
                infectedState[patient] = true;
                wasEverInfected[patient] = true;
            }

            //build and write heat map of infections
            std::string filename = std::string("heatMap-") + country + std::string("-step-") + std::to_string(step) + std::string(".dat");
            convertToHeatMap(infectedState, xcoords, ycoords, filename, resolution);

            step++;
            if (step == 5000) {
                runtime.stop();
                INFO("5000 steps took:", runtime.elapsedTag());
            }
        }
        INFO("Total infections: ", std::count(wasEverInfected.begin(), wasEverInfected.end(), true));
    }
}
} /* namespace NetworKit */
