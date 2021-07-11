// no-networkit-format
/*
 * HyperbolicGenerator.hpp
 *
 *  Created on: 20.05.2014
 *      Author: Moritz v. Looz
 */

#ifndef NETWORKIT_GENERATORS_HYPERBOLIC_GENERATOR_HPP_
#define NETWORKIT_GENERATORS_HYPERBOLIC_GENERATOR_HPP_

#include <cmath>
#include <vector>

#include <networkit/auxiliary/Timer.hpp>
#include <networkit/geometric/HyperbolicSpace.hpp>
#include <networkit/generators/StaticGraphGenerator.hpp>
#include <networkit/generators/quadtree/Quadtree.hpp>

namespace NetworKit {

/**
 * @ingroup generators
 */
class HyperbolicGenerator final : public StaticGraphGenerator {
    friend class DynamicHyperbolicGenerator;
public:

    /**
     * @param[in] n Number of nodes
     * @param[in] k Target average degree
     * @param[in] exp Target exponent of power-law distribution
     * @param[in] T Temperature
     */
    HyperbolicGenerator(count n=10000, double avgDegree=6, double exp=3, double T=0);

    /**
     * @param[in] angles Pointer to angles of node positions
     * @param[in] radii Pointer to radii of node positions
     * @param[in] r radius of poincare disk to place nodes in
     * @param[in] thresholdDistance Edges are added for nodes closer to each other than this threshold
     * @return Graph to be generated according to parameters
     */
    Graph generate(const vector<double> &angles, const vector<double> &radii, double R, double T=0);
    Graph generateCold(const vector<double> &angles, const vector<double> &radii, double R);

    /**
     * @return Graph to be generated according to parameters specified in constructor.
     */
    Graph generate() override;

    /**
     * Set the capacity of a quadtree leaf.
     *
     * @param capacity Tuning parameter, default value is 1000
     */
    void setLeafCapacity(count capacity) {
        this->capacity = capacity;
    }

    /**
     * When using a theoretically optimal split, the quadtree will be flatter, but running time usually longer.
     * @param split Whether to use the theoretically optimal split. Defaults to false
     */
    void setTheoreticalSplit(bool split) {
        this->theoreticalSplit = split;
    }

    void setBalance(double balance) {
        this->balance = balance;
    }

    vector<double> getElapsedMilliseconds() const {
        vector<double> result(threadtimers.size());
        for (index i = 0; i < result.size(); i++) {
            result[i] = threadtimers[i].elapsedMilliseconds();
        }
        return result;
    }

private:

    /**
     * Set tuning parameters to their default values
     */
    void initialize();

    Graph generate(count n, double R, double alpha, double T = 0);

    static vector<vector<double>> getBandAngles(const vector<vector<Point2DWithIndex<double>>> &bands) {
        vector<vector<double>> bandAngles(bands.size());
        #pragma omp parallel for
        for (omp_index i=0; i < static_cast<omp_index>(bands.size()); i++){
            const count currentBandSize = bands[i].size();
            bandAngles[i].resize(currentBandSize);
            for(index j=0; j < currentBandSize; j++) {
                bandAngles[i][j] = bands[i][j].getX();
            }
        }
        return bandAngles;
    }

    static vector<double> getBandRadii(int n, double R, double seriesRatio = 0.9) {
        /*
        * We assume band differences form a geometric series.
        * Thus, there is a constant ratio(r) between band length differences
        * i.e (c2-c1)/(c1-c0) = (c3-c2)/(c2-c1) = r
        */
        vector<double> bandRadius;
        bandRadius.push_back(0);
        double a = R*(1-seriesRatio)/(1-pow(seriesRatio, log(n)));
        const double logn = log(n);

        for (int i = 1; i < logn; i++){
            double c_i = a*((1-pow(seriesRatio, i))/(1-seriesRatio));
            bandRadius.push_back(c_i);
        }
        bandRadius.push_back(R);
        return bandRadius;
    }

    static std::tuple<double, double> getMinMaxTheta(double angle, double radius, double cLow, double thresholdDistance) {
      /*
          Calculates the angles that are enclosing the intersection of the
          hyperbolic disk that is around point v and the bands.
          Calculation is as follows:
          1. For the most inner band, return [0, 2pi]
          2. For other bands, consider the point P which lies on the tangent from origin to the disk of point v.
          Its radial coordinates would be(cHigh, point[1]+deltaTheta). We're looking for the deltaTheta.
          We know the distance from point v to P is R. Thus, we can solve the hyperbolic distance of (v, P)
          for deltaTheta. Then, thetaMax is simply point[1] + deltaTheta and thetaMin is point[1] - deltaTheta
      */

      //Most innerband is defined by cLow = 0
      double minTheta, maxTheta;
      if (cLow == 0)
      return std::make_tuple(0.0, 2* PI);

      double a = (cosh(radius)*cosh(cLow) - cosh(thresholdDistance))/(sinh(radius)*sinh(cLow));
      //handle floating point error
      if(a < -1)
        a = -1;
      else if(a > 1)
        a = 1;
      a = acos(a);
      maxTheta = angle + a;
      minTheta = angle - a;
      return std::make_tuple(minTheta, maxTheta);
    }

    static vector<Point2DWithIndex<double>> getPointsWithinAngles(double minTheta, double maxTheta, const vector<Point2DWithIndex<double>> &band, vector<double> &bandAngles){
        /**
        Returns the list of points, w, that lies within minTheta and maxTheta
        in the supplied band(That area is called as slab)
        */
        assert(band.size() == bandAngles.size());

        vector<Point2DWithIndex<double>> slab;

        std::vector<double>::iterator low;
        std::vector<double>::iterator high;

        if(minTheta == -2*PI)
            minTheta = 0;
        //Case 1: We do not have overlap 2pi, simply put all the points between min and max to the list
        if(maxTheta <= 2*PI && minTheta >= 0){
            low = std::lower_bound(bandAngles.begin(), bandAngles.end(), minTheta);
            high = std::upper_bound(bandAngles.begin(), bandAngles.end(), maxTheta);
            auto first = band.begin() + (low - bandAngles.begin());
            auto last = band.begin() + (high - bandAngles.begin());
            //Q: Does this operation increases the complexity ? It is linear in times of high - low
            //Does not increase the complexity, since we have to check these points anyway
            slab.insert(slab.end(), first, last);
        }
        //Case 2: We have 'forward' overlap at 2pi, that is maxTheta > 2pi
        else if (maxTheta > 2*PI){
            //1. Get points from minTheta to 2pi
            low = std::lower_bound(bandAngles.begin(), bandAngles.end(), minTheta);
            high = std::upper_bound(bandAngles.begin(), bandAngles.end(), 2*PI);
            auto first = band.begin() + (low - bandAngles.begin());
            auto last = band.begin() + (high - bandAngles.begin());
            slab.insert(slab.end(), first, last);

            //2. Get points from 0 to maxTheta%2pi
            low = std::lower_bound(bandAngles.begin(), bandAngles.end(), 0);
            maxTheta = fmod(maxTheta, (2*PI));
            high = std::upper_bound(bandAngles.begin(), bandAngles.end(), maxTheta);
            auto first2 = band.begin() + (low - bandAngles.begin());
            auto last2 = band.begin() + (high - bandAngles.begin());
            slab.insert(slab.end(), first2, last2);
        }
        //Case 3: We have 'backward' overlap at 2pi, that is minTheta < 0
        else if (minTheta < 0){
            //1. Get points from 2pi + minTheta to 2pi
            minTheta = (2*PI) + minTheta;
            low = std::lower_bound(bandAngles.begin(), bandAngles.end(), minTheta);
            high = std::upper_bound(bandAngles.begin(), bandAngles.end(), 2*PI);
            auto first = band.begin() + (low - bandAngles.begin());
            auto last = band.begin() + (high - bandAngles.begin());
            slab.insert(slab.end(), first, last);

            //2. Get points from 0 to maxTheta
            low = std::lower_bound(bandAngles.begin(), bandAngles.end(), 0);
            high = std::upper_bound(bandAngles.begin(), bandAngles.end(), maxTheta);
            auto first2 = band.begin() + (low - bandAngles.begin());
            auto last2 = band.begin() + (high - bandAngles.begin());
            slab.insert(slab.end(), first2, last2);
        }
        return slab;
    }

    /**
     * graph parameters
     */
    count nodeCount;
    double R;
    double alpha;
    double temperature;

    /**
     * tuning parameters
     */
    count capacity;
    bool theoreticalSplit;
    double balance = 0.5;
    static constexpr bool directSwap = false;

    /**
     * times
     */
    vector<Aux::Timer> threadtimers;
};
}
#endif // NETWORKIT_GENERATORS_HYPERBOLIC_GENERATOR_HPP_
