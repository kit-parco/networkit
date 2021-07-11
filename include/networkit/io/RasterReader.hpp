/*
 * RasterReader.hpp
 *
 *  Created on: 20.07.2015
 *      Author: Henning
 */

#ifndef NETWORKIT_IO_RASTER_READER_HPP_
#define NETWORKIT_IO_RASTER_READER_HPP_

#include <string>
#include <vector>

namespace NetworKit {

/**
 * Reader for NASA raster data of population.
 * Coordinate transformation into the unit square takes place.
 * Within each grid cell, points get a random location.
 */
class RasterReader final {
    double normalizationFactor;

public:
    /**
     * @param[in] Factor by which population per grid cell should be multiplied
     */
    RasterReader(double normalizationFactor);

    ~RasterReader() = default;

    /**
     * @param[in] path Path and name of file to be read.
     * @return Vectors of x- and y-coordinates of the points.
     */
    std::pair<std::vector<double>, std::vector<double>> read(const std::string &path);
};

} /* namespace NetworKit */
#endif // NETWORKIT_IO_RASTER_READER_HPP_
