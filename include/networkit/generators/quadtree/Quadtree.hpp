// no-networkit-format
/*
 * Quadtree.hpp
 *
 *  Created on: 21.05.2014
 *      Author: Moritz v. Looz
 */

#ifndef NETWORKIT_GENERATORS_QUADTREE_QUADTREE_HPP_
#define NETWORKIT_GENERATORS_QUADTREE_QUADTREE_HPP_

#include <cmath>
#include <memory>
#include <omp.h>
#include <vector>

#include <networkit/auxiliary/Parallel.hpp>
#include <networkit/generators/quadtree/QuadNode.hpp>
#include <networkit/geometric/HyperbolicSpace.hpp>

namespace NetworKit {

template <class T, bool poincare=true>
class Quadtree final {
    friend class QuadTreeGTest;
public:
    Quadtree() {
        root = QuadNode<T, poincare>();
        this->maxRadius = 1;
    }

    /**
     * @param maxR Radius of the managed area. Must be smaller than 1.
     * @param theoreticalSplit If true, split cells to get the same area in each child cell. Default is false
     * @param alpha dispersion Parameter of the point distribution. Only has an effect if theoretical split is true
     * @param capacity How many points can inhabit a leaf cell before it is split up?
     *
     */
    Quadtree(double maxR,bool theoreticalSplit=false, double alpha=1, count capacity=1000, double balance = 0.5) {
        root = QuadNode<T,poincare>(0, 0, 2*PI, maxR, capacity, theoreticalSplit,alpha,balance);
        this->maxRadius = maxR;
    }

    Quadtree(const vector<double> &angles, const vector<double> &radii, const vector<T> &content, double stretch, bool theoreticalSplit=false, double alpha=1, count capacity=1000, double balance = 0.5) {
        const count n = angles.size();
        assert(angles.size() == radii.size());
        assert(radii.size() == content.size());
        double R = stretch*HyperbolicSpace::hyperbolicAreaToRadius(n);
        double r = HyperbolicSpace::hyperbolicRadiusToEuclidean(R);
        root = QuadNode<T>(0, 0, 2*PI, r, capacity, theoreticalSplit,alpha,balance);
        maxRadius = r;
        for (index i = 0; i < n; i++) {
            assert(content[i] < n);
            root.addContent(content[i], angles[i], radii[i]);
        }
    }

    /**
     * @param newcomer content to be added at point x
     * @param angle angular coordinate of x
     * @param R radial coordinate of x
     */
    void addContent(T newcomer, double angle, double r) {
        root.addContent(newcomer, angle, r);
    }

    /**
     * @param newcomer content to be removed at point x
     * @param angle angular coordinate of x
     * @param R radial coordinate of x
     */
    bool removeContent(T toRemove, double angle, double r) {
        return root.removeContent(toRemove, angle, r);
    }

    /**
     * Get all elements, regardless of position
     *
     * @return vector<T> of elements
     */
    vector<T> getElements() const {
        return root.getElements();
    }

    void extractCoordinates(vector<double> &anglesContainer, vector<double> &radiiContainer) const {
        root.getCoordinates(anglesContainer, radiiContainer);
    }

    /**
     * Get elements whose hyperbolic distance to the query point is less than the hyperbolic distance
     *
     * @param circleCenter Cartesian coordinates of the query circle's center
     * @param hyperbolicRadius Radius of the query circle
     */
    vector<T> getElementsInHyperbolicCircle(Point2DWithIndex<double> circleCenter, double hyperbolicRadius) const {
        vector<T> circleDenizens;
        getElementsInHyperbolicCircle(circleCenter, hyperbolicRadius, circleDenizens);
        return circleDenizens;
    }

    void getElementsInHyperbolicCircle(const Point2DWithIndex<double> circleCenter, const double hyperbolicRadius, const bool suppressLeft, vector<T> &circleDenizens) const {
        assert(circleDenizens.empty());
        double cc_phi, cc_r;
        HyperbolicSpace::cartesianToPolar(circleCenter, cc_phi, cc_r);
        //Transform hyperbolic circle into Euclidean circle
        double minPhi, maxPhi, radius, r_e;
        HyperbolicSpace::getEuclideanCircle(cc_r, hyperbolicRadius, r_e, radius);
        Point2DWithIndex<double> center = HyperbolicSpace::polarToCartesian(cc_phi, r_e);
        double minR = r_e - radius;
        double maxR = r_e + radius;
        if (maxR > 1) maxR = 1;
        if (minR < 0) {
            maxR = std::max(abs(minR), maxR);
            minR = 0;
            minPhi = 0;
            maxPhi = 2*PI;
        } else {
            double spread = asin(radius / r_e);
            minPhi = cc_phi - spread;
            maxPhi = cc_phi + spread;
            /**
             * If the circle overlaps the 2\pi line, we have to make two separate calls and collect
             */
        }

        if (suppressLeft) minPhi = cc_phi;
        /**
         * get Elements in Euclidean circle
         */

        bool wraparound = false;
        root.getElementsInEuclideanCircle(center, radius, circleDenizens, minPhi, maxPhi, minR, maxR);
        if (minPhi < 0) {
            root.getElementsInEuclideanCircle(center, radius, circleDenizens, 2*PI+minPhi, 2*PI, minR, maxR);
            wraparound = true;
        }
        if (maxPhi > 2*PI) {
            root.getElementsInEuclideanCircle(center, radius, circleDenizens, 0, maxPhi - 2*PI, minR, maxR);
            wraparound = true;
        }

        for (T denizen : circleDenizens) {
            if (denizen >= size()) {
                DEBUG("Content ", denizen, " found in quadtree of size ", size(), ", as one of ", circleDenizens.size(), " neighbors.");
            }
        }

        //we have sort(deg(v)) here! This is not good, but does not make the asymptotical complexity of O(deg(v) log n) worse.
        if (wraparound) {
            Aux::Parallel::sort(circleDenizens.begin(), circleDenizens.end());
            auto newend = unique(circleDenizens.begin(), circleDenizens.end());
            count toRemove = circleDenizens.end() - newend;
            count remaining = newend - circleDenizens.begin();
            if (toRemove > 0) {
                DEBUG("Removing, ", toRemove, " duplicate entries, keeping ", remaining);
                circleDenizens.resize(remaining);
            }
        }

        for (T denizen : circleDenizens) {
            if (denizen >= size()) DEBUG("Content ", denizen, " found in quadtree of size ", size(), ", as one of ", circleDenizens.size(), " neighbors, after sorting");
        }
    }

    void getElementsInHyperbolicCircle(const Point2DWithIndex<double> circleCenter, const double hyperbolicRadius, vector<T> &circleDenizens) const {
        getElementsInHyperbolicCircle(circleCenter, hyperbolicRadius, false, circleDenizens);
    }

    count getElementsProbabilistically(Point2DWithIndex<double> euQuery, std::function<double(double)> prob, vector<T> &circleDenizens) {
        return root.getElementsProbabilistically(euQuery, prob, false, circleDenizens);
    }

    count getElementsProbabilistically(Point2DWithIndex<double> euQuery, std::function<double(double)> prob, bool suppressLeft, vector<T> &circleDenizens) {
        return root.getElementsProbabilistically(euQuery, prob, suppressLeft, circleDenizens);
    }

    void recount() {
        root.recount();
    }

    count size() const {
        return root.size();
    }

    count height() const {
        return root.height();
    }

    count countLeaves() const {
        return root.countLeaves();
    }

    index indexSubtree(index nextID) {
        return root.indexSubtree(nextID);
    }

    index getCellID(double phi, double r) const {
        return root.getCellID(phi, r);
    }

    double getMaxRadius() const {
        return maxRadius;
    }

    void reindex() {
        #pragma omp parallel
        {
            #pragma omp single nowait
            {
                root.reindex(0);
            }
        }
    }

    /**
     * trims the vectors used to hold the content in the leaf cells. Reduces memory usage, makes changes slower
     */
    void trim() {
        root.trim();
    }

private:
    QuadNode<T, poincare> root;
    double maxRadius;
};
}

#endif // NETWORKIT_GENERATORS_QUADTREE_QUADTREE_HPP_
