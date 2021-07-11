// no-networkit-format
/*
 * RasterReader.cpp
 *
 *  Created on: 20.07.2015
 *      Author: Henning
 */

#include <fstream>

#include <networkit/auxiliary/Log.hpp>
#include <networkit/auxiliary/NumberParsing.hpp>
#include <networkit/auxiliary/Random.hpp>
#include <networkit/io/RasterReader.hpp>

namespace NetworKit {

RasterReader::RasterReader(double normalizationFactor): normalizationFactor(normalizationFactor) {}

std::pair<std::vector<double>, std::vector<double>>
RasterReader::read(const std::string& path) {
    DEBUG("start reading raster file...");

    std::ifstream file(path);
    Aux::enforceOpened(file);
    std::string line;

    DEBUG("file opened: ", path);

    // *** read header: 6 lines, ignore 3-6
    double val = 0.0;

    // ** read number of columns
    std::getline(file, line);
    auto it = line.begin() + 14; // ignore the first 14 characters
    auto end = line.end();
    std::tie(val, it) = Aux::Parsing::strTo<double>(it, end);
    uint64_t ncols = (uint64_t) val;
    DEBUG("line: ", line, "; nrows: ", ncols);

    // ** read number of rows
    std::getline(file, line);
    it = line.begin() + 14; // ignore the first 14 characters
    end = line.end();
    std::tie(val, it) = Aux::Parsing::strTo<double>(it, end);
    uint64_t nrows = (uint64_t) val;
    DEBUG("line: ", line, "; ncols: ", nrows);

    // ** ignore remaining header info
    for (uint64_t i = 3; i <= 6; ++i) {
        std::getline(file, line);
    }


    // data for main loop
    double xdelta = 1.0;// / (double) ncols; // width of each column
    double ydelta = 1.0;// / (double) nrows; // height of each row
    double xlb = 0.0;                     // lower bound in x-direction
    double ylb = 1.0 - ydelta;            // lower bound in y-direction
    double xub = xdelta;                  // upper bound in x-direction
    double yub = 1.0;                     // upper bound in y-direction
    double validate = 0.0;                // to validate the sum of the values

    DEBUG("xdelta: ", xdelta, ", ydelta: ", ydelta);

    std::vector<double> xcoords;          // stores the x-coordinates of generated points
    std::vector<double> ycoords;          // stores the y-coordinates of generated points

    for (uint64_t row = 0; row < nrows; ++row) {
        std::getline(file, line);
        it = line.begin();
        end = line.end();

        for (uint64_t col = 0; col < ncols; ++col) {
            // read next number
            std::tie(val, it) = Aux::Parsing::strTo<double>(it, end);
            validate += val;

            // divide by "normalizer" and round down
            val *= normalizationFactor;
            uint64_t numPointsInCell = (uint64_t) std::round(val);

            for (uint64_t i = 0; i < numPointsInCell; ++i) {
                // insert random coordinate into coordinate arrays
                xcoords.push_back(Aux::Random::real(xlb, xub));
                ycoords.push_back(Aux::Random::real(ylb, yub));
            }

            // adjust x-coordinate
            xlb = xub;
            xub += xdelta;
        }

        // adjust coordinates: reset x, decrease y
        xlb = 0.0;
        xub = xdelta;
        yub = ylb;
        ylb -= ydelta;
    }
    DEBUG("Created ", xcoords.size(), " 2D points, validation: ", validate, " vs ", xcoords.size() / normalizationFactor);
    return std::make_pair(xcoords, ycoords);
}

} /* namespace NetworKit */
