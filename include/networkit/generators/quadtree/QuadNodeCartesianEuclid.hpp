// no-networkit-format
/*
 * QuadNodePolarEuclid.hpp
 *
 *  Created on: 21.05.2014
 *      Author: Moritz v. Looz
 *
 *  Note: This is similar enough to QuadNode.hpp that one could merge these two classes.
 */

#ifndef NETWORKIT_GENERATORS_QUADTREE_QUAD_NODE_CARTESIAN_EUCLID_HPP_
#define NETWORKIT_GENERATORS_QUADTREE_QUAD_NODE_CARTESIAN_EUCLID_HPP_

#include <algorithm>
#include <cassert>
#include <functional>
#include <vector>

#include <networkit/auxiliary/Log.hpp>
#include <networkit/geometric/HyperbolicSpace.hpp>

using std::vector;
using std::min;
using std::max;
using std::cos;

namespace NetworKit {

template <class T>
class QuadNodeCartesianEuclid final {
    friend class QuadTreeGTest;

    Point<double> minPoint;
    Point<double> maxPoint;
    count dimension;
    unsigned capacity;
    static const unsigned coarsenLimit = 4;
    count subTreeSize;
    std::vector<T> content;
    std::vector<Point<double> > positions;
    bool isLeaf;
    bool splitTheoretical;
    index ID;
    double lowerBoundR;

public:
    std::vector<QuadNodeCartesianEuclid> children;

    /**
     * Construct a QuadNode for polar coordinates.
     *
     *
     * @param lower Minimal coordinates of region
     * @param upper Maximal coordinates of region (excluded)
     * @param capacity Number of points a leaf cell can store before splitting
     * @param splitTheoretical Whether to split in a theoretically optimal way or in a way to decrease measured running times
     *
     */
    QuadNodeCartesianEuclid(const Point<double> &lower = Point<double>(0.0, 0.0),
                            const Point<double> &upper = Point<double>(1.0, 1.0),
                            unsigned capacity = 1000, bool splitTheoretical = false) {
        this->minPoint = lower;
        this->maxPoint = upper;
        this->dimension = minPoint.getDimensions();
        assert(maxPoint.getDimensions() == dimension);
        this->capacity = capacity;
        this->splitTheoretical = splitTheoretical;
        this->ID = 0;
        isLeaf = true;
        subTreeSize = 0;
    }

    void split() {
        assert(isLeaf);
        assert(children.size() == 0);
        vector<double> middle(dimension);
        if (splitTheoretical) {
            //Euclidean space is distributed equally
            for (index d = 0; d < dimension; d++) {
                middle[d] = (minPoint[d] + maxPoint[d]) / 2;
            }
        } else {
            //median of points
            const count numPoints = positions.size();
            assert(numPoints > 0);//otherwise, why split?
            vector<vector<double> > sorted(dimension);
            for (index d = 0; d < dimension; d++) {
                sorted[d].resize(numPoints);
                for (index i = 0; i < numPoints; i++) {
                    sorted[d][i] = positions[i][d];
                }
                std::sort(sorted[d].begin(), sorted[d].end());
                middle[d] = sorted[d][numPoints/2];//this will crash if no points are there!
                assert(middle[d] <= maxPoint[d]);
                assert(middle[d] >= minPoint[d]);
            }
        }
        count childCount = pow(2,dimension);
        for (index i = 0; i < childCount; i++) {
            vector<double> lowerValues(dimension);
            vector<double> upperValues(dimension);
            index bitCopy = i;
            for (index d = 0; d < dimension; d++) {
                if (bitCopy & 1) {
                    lowerValues[d] = middle[d];
                    upperValues[d] = maxPoint[d];
                } else {
                    lowerValues[d] = minPoint[d];
                    upperValues[d] = middle[d];
                }
                bitCopy = bitCopy >> 1;
            }
            QuadNodeCartesianEuclid child(Point<double>(lowerValues), Point<double>(upperValues), capacity, splitTheoretical);
            assert(child.isLeaf);
            children.push_back(child);
        }
        isLeaf = false;
    }

    /**
     * Add a point at polar coordinates (angle, R) with content input. May split node if capacity is full
     *
     * @param input arbitrary content, in our case an index
     * @param angle angular coordinate of point, between 0 and 2 pi.
     * @param R radial coordinate of point, between 0 and 1.
     */
    void addContent(T input, Point<double> pos) {
        assert(content.size() == positions.size());
        assert(this->responsible(pos));
        if (isLeaf) {
            if (content.size() + 1 < capacity) {
                content.push_back(input);
                positions.push_back(pos);
            } else {
                split();

                for (index i = 0; i < content.size(); i++) {
                    this->addContent(content[i], positions[i]);
                }
                assert(subTreeSize == content.size());//we have added everything twice
                subTreeSize = content.size();
                content.clear();
                positions.clear();
                this->addContent(input, pos);
            }
        }
        else {
            assert(children.size() > 0);
            bool foundResponsibleChild = false;
            for (index i = 0; i < children.size(); i++) {
                if (children[i].responsible(pos)) {
                    foundResponsibleChild = true;
                    children[i].addContent(input, pos);
                    break;
                }
            }
            assert(foundResponsibleChild);
            (void)foundResponsibleChild;
            subTreeSize++;
        }
    }

    /**
     * Remove content at coordinate pos. May cause coarsening of the quadtree
     *
     * @param input Content to be removed
     * @param pos Coordinate of content
     *
     * @return True if content was found and removed, false otherwise
     */
    bool removeContent(T input, Point<double> pos) {
        if (!responsible(pos)) return false;
        if (isLeaf) {
            index i = 0;
            for (; i < content.size(); i++) {
                if (content[i] == input) break;
            }
            if (i < content.size()) {
                assert(positions[i].distance(pos) == 0);
                //remove element
                content.erase(content.begin()+i);
                positions.erase(positions.begin()+i);
                return true;
            } else {
                return false;
            }
        }
        else {
            bool removed = false;
            bool allLeaves = true;
            assert(children.size() > 0);
            for (index i = 0; i < children.size(); i++) {
                if (!children[i].isLeaf) allLeaves = false;
                if (children[i].removeContent(input, pos)) {
                    assert(!removed);
                    removed = true;
                }
            }
            if (removed) subTreeSize--;
            //coarsen?
            if (removed && allLeaves && size() < coarsenLimit) {
                //coarsen!!
                //why not assert empty containers and then insert directly?
                vector<T> allContent;
                vector<Point<double> > allPositions;
                for (index i = 0; i < children.size(); i++) {
                    allContent.insert(allContent.end(), children[i].content.begin(), children[i].content.end());
                    allPositions.insert(allPositions.end(), children[i].positions.begin(), children[i].positions.end());
                }
                assert(allContent.size() == allPositions.size());
                children.clear();
                content.swap(allContent);
                positions.swap(allPositions);
                isLeaf = true;
            }

            return removed;
        }
    }


    /**
     * Check whether the region managed by this node lies outside of an Euclidean circle.
     *
     * @param query Center of the Euclidean query circle, given in Cartesian coordinates
     * @param radius Radius of the Euclidean query circle
     *
     * @return True if the region managed by this node lies completely outside of the circle
     */
    bool outOfReach(const Point<double> &query, double radius) const {
        return EuclideanDistances(query).first > radius;
    }

    /**
     * @param query Position of the query point
     */
    std::pair<double, double> EuclideanDistances(const Point<double> &query) const {
        /**
         * If the query point is not within the quadnode, the distance minimum is on the border.
         * Need to check whether extremum is between corners.
         */
        double maxDistance = 0;
        double minDistance = std::numeric_limits<double>::max();

        if (responsible(query)) minDistance = 0;

        auto updateMinMax = [&minDistance, &maxDistance, &query](const Point<double> &pos) -> void {
            double extremalValue = pos.distance(query);
            maxDistance = std::max(extremalValue, maxDistance);
            minDistance = std::min(minDistance, extremalValue);
        };

        vector<double> closestValues(dimension);
        vector<double> farthestValues(dimension);

        for (index d = 0; d < dimension; d++) {
            if (std::abs(query[d] - minPoint.at(d)) < std::abs(query[d] - maxPoint.at(d))) {
                closestValues[d] = minPoint.at(d);
                farthestValues[d] = maxPoint.at(d);
            } else {
                farthestValues[d] = minPoint.at(d);
                closestValues[d] = maxPoint.at(d);
            }
            if (query[d] >= minPoint.at(d) && query[d] <= maxPoint.at(d)) {
                closestValues[d] = query[d];
            }
        }
        updateMinMax(Point<double>(closestValues));
        updateMinMax(Point<double>(farthestValues));

        assert(minDistance < query.length() + maxPoint.length());
        assert(minDistance < maxDistance);
        return {minDistance, maxDistance};
    }


    /**
     * Does the point at (angle, r) fall inside the region managed by this QuadNode?
     *
     * @param angle Angular coordinate of input point
     * @param r Radial coordinate of input points
     *
     * @return True if input point lies within the region of this QuadNode
     */
    bool responsible(Point<double> pos) const {
        for (index d = 0; d < dimension; d++) {
            if (pos[d] < minPoint.at(d) || pos[d] >= maxPoint.at(d)) return false;
        }
        return true;
    }

    /**
     * Get all Elements in this QuadNode or a descendant of it
     *
     * @return vector of content type T
     */
    std::vector<T> getElements() const {
        if (isLeaf) {
            return content;
        } else {
            assert(content.size() == 0);
            assert(positions.size() == 0);
            vector<T> result;
            for (index i = 0; i < children.size(); i++) {
                std::vector<T> subresult = children[i].getElements();
                result.insert(result.end(), subresult.begin(), subresult.end());
            }
            return result;
        }
    }

    void getCoordinates(vector<Point<double> > &pointContainer) const {
        if (isLeaf) {
            pointContainer.insert(pointContainer.end(), positions.begin(), positions.end());
        }
        else {
            assert(content.size() == 0);
            assert(positions.size() == 0);
            for (index i = 0; i < children.size(); i++) {
                children[i].getCoordinates(pointContainer);
            }
        }
    }

    /**
     * Main query method, get points lying in a Euclidean circle around the center point.
     * Optional limits can be given to get a different result or to reduce unnecessary comparisons
     *
     * Elements are pushed onto a vector which is a required argument. This is done to reduce copying.
     * (Maybe not necessary due to copy elisison)
     *
     * Safe to call in parallel.
     *
     * @param center Center of the query circle
     * @param radius Radius of the query circle
     * @param result Reference to the vector where the results will be stored
     * @param minAngle Optional value for the minimum angular coordinate of the query region
     * @param maxAngle Optional value for the maximum angular coordinate of the query region
     * @param lowR Optional value for the minimum radial coordinate of the query region
     * @param highR Optional value for the maximum radial coordinate of the query region
     */
    void getElementsInEuclideanCircle(Point<double> center, double radius, vector<T> &result) const {
        if (outOfReach(center, radius)) {
            return;
        }

        if (isLeaf) {
            const double rsq = radius*radius;
            const count cSize = content.size();

            for (index i=0; i < cSize; i++) {
                if (positions[i].squaredDistance(center) < rsq) {
                    result.push_back(content[i]);
                }
            }
        }  else {
            for (index i = 0; i < children.size(); i++) {
                children[i].getElementsInEuclideanCircle(center, radius, result);
            }
        }
    }

    count getElementsProbabilistically(Point<double> euQuery, std::function<double(double)> prob, vector<T> &result) const {
        TRACE("Getting Euclidean distances");
        auto distancePair = EuclideanDistances(euQuery);
        double probUB = prob(distancePair.first);
        double probLB = prob(distancePair.second);
        assert(probLB <= probUB);
        if (probUB > 0.5) probUB = 1;
        if (probUB == 0) return 0;
        //TODO: return whole if probLB == 1
        double probdenom = std::log(1-probUB);
        if (probdenom == 0) return 0;//there is a very small probability, but we cannot process it.
        TRACE("probUB: ", probUB, ", probdenom: ", probdenom);

        count expectedNeighbours = probUB*size();
        count candidatesTested = 0;
        count incomingNeighbours = result.size();
        count ownsize = size();


        if (isLeaf) {
            const count lsize = content.size();
            TRACE("Leaf of size ", lsize);
            for (index i = 0; i < lsize; i++) {
                //jump!
                if (probUB < 1) {
                    double random = Aux::Random::real();
                    double delta = std::log(random) / probdenom;
                    assert(delta >= 0);
                    i += delta;
                    if (i >= lsize) break;
                    TRACE("Jumped with delta ", delta, " arrived at ", i);
                }

                //see where we've arrived
                candidatesTested++;
                double distance = positions[i].distance(euQuery);
                assert(distance >= distancePair.first);//TODO: These should not fail!
                assert(distance <= distancePair.second);
                double q = prob(distance);
                q = q / probUB; //since the candidate was selected by the jumping process, we have to adjust the probabilities
                assert(q <= 1);

                //accept?
                double acc = Aux::Random::real();
                if (acc < q) {
                    TRACE("Accepted node ", i, " with probability ", q, ".");
                    result.push_back(content[i]);
                }
            }
        }  else {
            if (expectedNeighbours < 4 || probUB < 1/1000) {//select candidates directly instead of calling recursively
                TRACE("probUB = ", probUB,  ", switching to direct candidate selection.");
                assert(probUB < 1);
                const count stsize = size();
                for (index i = 0; i < stsize; i++) {
                    double delta = std::log(Aux::Random::real()) / probdenom;
                    assert(delta >= 0);
                    i += delta;
                    TRACE("Jumped with delta ", delta, " arrived at ", i, ". Calling maybeGetKthElement.");
                    if (i < size()) maybeGetKthElement(probUB, euQuery, prob, i, result);//this could be optimized. As of now, the offset is subtracted separately for each point
                    else break;
                    candidatesTested++;
                }
            } else {//carry on as normal
                for (index i = 0; i < children.size(); i++) {
                    TRACE("Recursively calling child ", i);
                    candidatesTested += children[i].getElementsProbabilistically(euQuery, prob, result);
                }
            }
        }
        count finalNeighbours = result.size();
        if (probLB == 1){
            assert(finalNeighbours == incomingNeighbours + ownsize);
            (void)finalNeighbours;
            (void)incomingNeighbours;
            (void)ownsize;
        }
        return candidatesTested;
    }


    void maybeGetKthElement(double upperBound, Point<double> euQuery, std::function<double(double)> prob, index k, vector<T> &circleDenizens) const {
        TRACE("Maybe get element ", k, " with upper Bound ", upperBound);
        assert(k < size());
        if (isLeaf) {
            double acceptance = prob(euQuery.distance(positions[k]))/upperBound;
            TRACE("Is leaf, accept with ", acceptance);
            if (Aux::Random::real() < acceptance) circleDenizens.push_back(content[k]);
        } else {
            TRACE("Call recursively.");
            index offset = 0;
            for (index i = 0; i < children.size(); i++) {
                count childsize = children[i].size();
                if (k - offset < childsize) {
                    children[i].maybeGetKthElement(upperBound, euQuery, prob, k - offset, circleDenizens);
                    break;
                }
                offset += childsize;
            }
        }
    }

    /**
     * Shrink all vectors in this subtree to fit the content.
     * Call after quadtree construction is complete, causes better memory usage and cache efficiency
     */
    void trim() {
        content.shrink_to_fit();
        positions.shrink_to_fit();
        if (!isLeaf) {
            for (index i = 0; i < children.size(); i++) {
                children[i].trim();
            }
        }
    }

    /**
     * Number of points lying in the region managed by this QuadNode
     */
    count size() const {
        return isLeaf ? content.size() : subTreeSize;
    }

    void recount() {
        subTreeSize = 0;
        for (index i = 0; i < children.size(); i++) {
            children[i].recount();
            subTreeSize += children[i].size();
        }
    }

    /**
     * Height of subtree hanging from this QuadNode
     */
    count height() const {
        count result = 1;//if leaf node, the children loop will not execute
        for (auto child : children) result = std::max(result, child.height()+1);
        return result;
    }

    /**
     * Leaf cells in the subtree hanging from this QuadNode
     */
    count countLeaves() const {
        if (isLeaf) return 1;
        count result = 0;
        for (index i = 0; i < children.size(); i++) {
            result += children[i].countLeaves();
        }
        return result;
    }

    index getID() const {
        return ID;
    }

    index indexSubtree(index nextID) {
        index result = nextID;
        assert(children.size() == pow(2,dimension) || children.size() == 0);
        for (int i = 0; i < children.size(); i++) {
            result = children[i].indexSubtree(result);
        }
        this->ID = result;
        return result+1;
    }

    index getCellID(Point<double> pos) const {
        if (!responsible(pos)) return none;
        if (isLeaf) return getID();
        else {
            for (int i = 0; i < children.size(); i++) {
                index childresult = children[i].getCellID(pos);
                if (childresult != none) return childresult;
            }
            throw std::runtime_error("No responsible child node found even though this node is responsible.");
        }
    }

    index getMaxIDInSubtree() const {
        if (isLeaf) return getID();
        else {
            index result = -1;
            for (int i = 0; i < children.size(); i++) {
                result = std::max(children[i].getMaxIDInSubtree(), result);
            }
            return std::max(result, getID());
        }
    }

    count reindex(count offset) {
        if (isLeaf)
        {
#ifndef NETWORKIT_OMP2
            #pragma omp task
#endif // NETWORKIT_OMP2
            {
                index p = offset;
                std::generate(content.begin(), content.end(), [&p](){return p++;});
            }
            offset += size();
        } else {
            for (int i = 0; i < children.size(); i++) {
                offset = children[i].reindex(offset);
            }
        }
        return offset;
    }
};
}

#endif // NETWORKIT_GENERATORS_QUADTREE_QUAD_NODE_CARTESIAN_EUCLID_HPP_
