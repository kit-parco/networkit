// no-networkit-format
/*
 * HyperbolicSpace.hpp
 *
 *  Created on: 20.05.2014
 *      Author: Moritz v. Looz
 */

#ifndef NETWORKIT_GEOMETRIC_HYPERBOLIC_SPACE_HPP_
#define NETWORKIT_GEOMETRIC_HYPERBOLIC_SPACE_HPP_

#include <cmath>
#include <map>
#include <vector>

#include <networkit/auxiliary/Random.hpp>
#include <networkit/geometric/Point2DWithIndex.hpp>
#include <networkit/viz/Point.hpp>

using std::vector;
using std::abs;

namespace NetworKit {

class HyperbolicSpace final {
public:
    HyperbolicSpace() = default;
    ~HyperbolicSpace() = default;
    /**
     * @param angles empty vector to hold angular coordinates of generated points
     * @param radii empty vector to hold radial coordinates of generated points
     * @param R radius of the hyperbolic disk
     * @param alpha dispersion parameter for the node positions
     *
     * Fill preallocated vectors with randomly sampled points in native coordinates
     */
    static void fillPoints(vector<double> &angles, vector<double> &radii, double R, double alpha);

    /**
     * @param angles empty vector to hold angular coordinates of generated points
     * @param radii empty vector to hold radial coordinates of generated points
     * @param stretch multiplier for the radius of the hyperbolic disk
     * @param alpha dispersion parameter for the node positions
     *
     * Fill preallocated vectors with randomly sampled points in native coordinates
     */
    static void fillPoints(vector<double> &angles, vector<double> &radii, double minPhi, double maxPhi, double minR, double maxR, double alpha);

    /**
     * @param firstangle angular coordinate of the first point
     * @param firstR radial coordinate of the first point
     * @param secondangle angular coordinate of the second point
     * @param secondR radial coordinate of the second point
     *
     * @return distance between two points in the poincare metric
     */
    static double poincareMetric(double firstangle, double firstR, double secondangle, double secondR);

    /**
     * @param firstangle angular coordinate of the first point
     * @param firstR radial coordinate of the first point
     * @param secondangle angular coordinate of the second point
     * @param secondR radial coordinate of the second point
     *
     * @return distance between two points in native hyperbolic representation
     */
    static double nativeDistance(double firstangle, double firstR, double secondangle, double secondR);

    /**
     * @param a first point in cartesian coordinates
     * @param b second point in cartesian coordinates
     *
     * @return distance between a and b in the poincare metric
     */
    static double poincareMetric(Point2DWithIndex<double> a, Point2DWithIndex<double> b);



    /**
     * @param phi angular coordinate of point
     * @param r radial coordinate of point
     *
     * @return cartesian coordinates represented by phi and r
     */
    static Point2DWithIndex<double> polarToCartesian(double phi, double r);

    /**
     * Convenience function for visualizations which expect coordinates as map<index,Point<float> >
     */
    static std::map<index, Point<float> > polarToCartesian(const vector<double> &angles, const vector<double> &radii);

    /**
     * @param a cartesian coordinates
     * @param phi empty double value to receive angular coordinate
     * @param r empty double value to receive radial coordinate
     */
    static void cartesianToPolar(Point2DWithIndex<double> a, double &phi, double &r);

    /**
     * Converts a hyperbolic circle to a Euclidean circle
     *
     * @param hyperbolicCenter center of the hyperbolic circle, given in cartesian coordinates within the poincare disk
     * @param hyperbolicRadius radius of the hyperbolic circle
     * @param euclideanCenter point to receive the center of the Euclidean circle, given in cartesian coordinates
     * @param euclidenRadius double to receive the radius of the Euclidean circle
     */
    static void getEuclideanCircle(Point2DWithIndex<double> hyperbolicCenter, double hyperbolicRadius, Point2DWithIndex<double> &euclideanCenter, double &euclideanRadius);
    static void getEuclideanCircle(double r_h, double hyperbolicRadius, double &radialCoordOfEuclideanCenter, double &euclideanRadius);

    /**
     * Project radial coordinates of the hyperbolic plane into the Poincare disk model
     *
     * @param hyperbolicRadius radial coordinate of a point in the native hyperbolic disc
     * @return radial coordinate in the Poincare model
     */
    static double hyperbolicRadiusToEuclidean(double hyperbolicRadius);

    /**
     * Project radial coordinates of the Poincare model into the hyperbolic plane
     *
     * @param EuclideanRadius Radial coordinate of a point in the Poincare model
     * @return radial coordinate in the hyperbolic plane
     */
    static double EuclideanRadiusToHyperbolic(double EuclideanRadius);

    /**
     * @param area The area of the hyperbolic circle
     * @param return Radius of a hyperbolic circle with the given area
     */
    static inline double hyperbolicAreaToRadius(double area) {
        return acosh(area/(2*PI)+1);
    }

    static inline double radiusToHyperbolicArea(double radius) {
            return  2*PI*(cosh(radius)-1);
    }

    static double getExpectedDegree(double n, double alpha, double R) {
        double gamma = 2*alpha+1;
        double xi = (gamma-1)/(gamma-2);
        double firstSumTerm = exp(-R/2);
        double secondSumTerm = exp(-alpha*R)*(alpha*(R/2)*((PI/4)*pow((1/alpha),2)-(PI-1)*(1/alpha)+(PI-2))-1);
        double expectedDegree = (2/PI)*xi*xi*n*(firstSumTerm + secondSumTerm);
        return expectedDegree;
    }

    static double searchTargetRadiusForColdGraphs(double n, double k, double alpha, double epsilon) {
        double gamma = 2*alpha+1;
        double xiInv = ((gamma-2)/(gamma-1));
        double v = k * (PI/2)*xiInv*xiInv;
        double currentR = 2*log(n / v);
        double lowerBound = currentR/2;
        double upperBound = currentR*2;
        assert(getExpectedDegree(n, alpha, lowerBound) > k);
        assert(getExpectedDegree(n, alpha, upperBound) < k);
        do {
            currentR = (lowerBound + upperBound)/2;
            double currentK = getExpectedDegree(n, alpha, currentR);
            if (currentK < k) {
                upperBound = currentR;
            } else {
                lowerBound = currentR;
            }
        } while (abs(getExpectedDegree(n, alpha, currentR) - k) > epsilon );
        return currentR;
    }

    static double getTargetRadius(double n, double m, double alpha=1, double T=0, double epsilon = 0.01) {
        double result;
        double plexp = 2*alpha+1;
        double targetAvgDegree = (m/n)*2;
        double xiInv = ((plexp-2)/(plexp-1));
        if (T == 0) {
            result = searchTargetRadiusForColdGraphs(n, targetAvgDegree, alpha, epsilon);
        } else {
            double beta = 1/T;
            if (T < 1){//cold regime
                double Iinv = ((beta/PI)*sin(PI/beta));
                double v = (targetAvgDegree*Iinv)*(PI/2)*xiInv*xiInv;
                result = 2*log(n / v);
            } else {//hot regime
                double v = targetAvgDegree*(1-beta)*pow((PI/2), beta)*xiInv*xiInv;
                result = 2*log(n/v)/beta;
            }
        }
        return result;
    }

    static inline double effectiveAreaInCell(double minPhi, double maxPhi, double minR, double maxR, double alpha) {
        double deltaPhi = maxPhi - minPhi;
        assert(deltaPhi >= 0);
        assert(deltaPhi <= 2*PI);
        double ringArea = radiusToHyperbolicArea(alpha*EuclideanRadiusToHyperbolic(maxR)) - radiusToHyperbolicArea(alpha*EuclideanRadiusToHyperbolic(minR));
        return ringArea*(deltaPhi/(2*PI));
    }

    /**
     * @param r_c radial coordinate of the circle center
     * @param d_c radius of the Euclidean circle
     * @param R radius of the hyperbolic base disk
     */
    static double hyperbolicSpaceInEuclideanCircle(double r_c, double d_c, double R);

    static double maxRinSlice(double minPhi, double maxPhi, double phi_c, double r_c, double euRadius);
};
}

#endif // NETWORKIT_GEOMETRIC_HYPERBOLIC_SPACE_HPP_
