// no-networkit-format
/*
 * QuadTreeGTest.cpp
 *
 *  Created on: 28.05.2014
 *      Author: Moritz v. Looz (moritz.looz-corswarem@kit.edu)
 */

#include <stack>
#include <cmath>
#include <algorithm>

#include <gtest/gtest.h>
#include <cmath>
#include <vector>
#include <networkit/generators/quadtree/Quadtree.hpp>

#include <networkit/auxiliary/Random.hpp>
#include <networkit/auxiliary/Log.hpp>
#include <networkit/geometric/HyperbolicSpace.hpp>

#include <networkit/generators/quadtree/QuadtreeCartesianEuclid.hpp>
#include <networkit/generators/quadtree/QuadtreePolarEuclid.hpp>

#include <tlx/unused.hpp>

namespace NetworKit {

class QuadTreeGTest: public testing::Test {
protected:
    template <class T>
    QuadNode<T> getRoot(Quadtree<T> &tree) {
        return tree.root;
    }

    template <class T>
    vector<QuadNode<T> > getChildren(QuadNode<T> &node) {
        return node.children;
    }
};


/**
 * Test whether the elements returned by a quadtree range query are indeed those whose hyperbolic distance to the query point is below a threshold
 */
TEST_F(QuadTreeGTest, testQuadTreeHyperbolicCircle) {
    count n = 1000;
    double R = 1;
    vector<double> angles(n);
    vector<double> radii(n);
    HyperbolicSpace::fillPoints(angles, radii, 1, 1);
    double max = 0;
    for (index i = 0; i < n; i++) {
        if (radii[i] > max) {
            max = radii[i];
        }
    }
    Quadtree<index> quad(max+(1-max)/4);

    for (index i = 0; i < n; i++) {
        EXPECT_GE(angles[i], 0);
        EXPECT_LT(angles[i], 2*PI);
        EXPECT_GE(radii[i], 0);
        EXPECT_LT(radii[i], R);
        TRACE("Added (", angles[i], ",", radii[i], ")");
        quad.addContent(i, angles[i], radii[i]);
    }
    vector<index> all = quad.getElements();
    EXPECT_EQ(all.size(), n);
    EXPECT_EQ(quad.size(), n);
    for (index testindex = 0; testindex < 100; testindex++) {
        index comparison = Aux::Random::integer(n-1);
        DEBUG("Using ", comparison, " at (", angles[comparison], ",", radii[comparison], ") as query point");

        vector<index> closeToOne = quad.getElementsInHyperbolicCircle(HyperbolicSpace::polarToCartesian(angles[comparison], radii[comparison]), R);
        EXPECT_LE(closeToOne.size(), n);

        for (index i = 0; i < closeToOne.size(); i++) {
            //no corrupt indices
            ASSERT_LE(closeToOne[i], n);

            //close results should actually be close
            EXPECT_LE(HyperbolicSpace::poincareMetric(angles[comparison], radii[comparison], angles[closeToOne[i]], radii[closeToOne[i]]), R);
            for (index j = 0; j < i; j++) {
                /**
                 * results are unique
                 */
                EXPECT_NE(closeToOne[i], closeToOne[j]);
            }
        }
        count notfound = 0;
        count didfind = 0;
        for (index i = 0; i < n; i++) {
            if (HyperbolicSpace::poincareMetric(angles[comparison], radii[comparison], angles[i], radii[i]) < R) {
                bool found = false;
                QuadNode<index> responsibleNode = getRoot(quad).getAppropriateLeaf(angles[i], radii[i]);

                for (index j = 0; j < closeToOne.size(); j++) {
                    if (closeToOne[j] == i) {
                        found = true;
                        break;
                    }
                }
                EXPECT_TRUE(found) << "dist(" << i << "," << comparison << ") = "
                        << HyperbolicSpace::poincareMetric(angles[comparison], radii[comparison], angles[i], radii[i]) << " < " << R;
                if (!found) {
                    notfound++;
                    DEBUG("angle: ", angles[i], ", radius: ", radii[i], ", leftAngle: ", responsibleNode.getLeftAngle(),
                            ", rightAngle: ", responsibleNode.getRightAngle(), ", minR: ", responsibleNode.getMinR(), ", maxR:", responsibleNode.getMaxR());
                    //DEBUG("euclidean Distance from circle center: ", center.distance(HyperbolicSpace::polarToCartesian(angles[i], radii[i])));
                    DEBUG("dist(", comparison, ", leftMin)=", HyperbolicSpace::poincareMetric(angles[comparison], radii[comparison], responsibleNode.getLeftAngle(), responsibleNode.getMinR()));
                    DEBUG("dist(", comparison, ", leftMax)=", HyperbolicSpace::poincareMetric(angles[comparison], radii[comparison], responsibleNode.getLeftAngle(), responsibleNode.getMaxR()));
                    DEBUG("dist(", comparison, ", rightMin)=", HyperbolicSpace::poincareMetric(angles[comparison], radii[comparison], responsibleNode.getRightAngle(), responsibleNode.getMinR()));
                    DEBUG("dist(", comparison, ", rightMax)=", HyperbolicSpace::poincareMetric(angles[comparison], radii[comparison], responsibleNode.getRightAngle(), responsibleNode.getMaxR()));
                }
                else {
                    didfind++;
                }
            }
        }
        if (notfound > 0) {
            DEBUG("Found only ", didfind, " of ", didfind + notfound, " neighbours");
        }
    }
}

/**
 * Gradually increase the distance threshold and check whether the number of neighbours increases monotonically. Necessary foundation for the dynamic hyperbolic generator.
 */
TEST_F(QuadTreeGTest, testQuadTreeThresholdGrowth) {
    count n = 100;
    double R = HyperbolicSpace::hyperbolicAreaToRadius(n);
    double r = HyperbolicSpace::hyperbolicRadiusToEuclidean(R);
    vector<double> angles(n);
    vector<double> radii(n);
    vector<double> indices(n);
    HyperbolicSpace::fillPoints(angles, radii, R, 1);
    double max = 0;
    for (index i = 0; i < n; i++) {
        indices[i] = i;
        radii[i] = HyperbolicSpace::hyperbolicRadiusToEuclidean(radii[i]);
        if (radii[i] == r) radii[i] = std::nextafter(radii[i], 0);
        if (radii[i] > max) {
            max = radii[i];
            assert(radii[i] < r);
        }
    }

    Quadtree<index> quad(r);

    for (index i = 0; i < n; i++) {
        EXPECT_GE(angles[i], 0);
        EXPECT_LT(angles[i], 2*PI);
        EXPECT_GE(radii[i], 0);
        EXPECT_LT(radii[i], R);
        TRACE("Added (", angles[i], ",", radii[i], ")");
        quad.addContent(i, angles[i], radii[i]);
    }

    for (index testindex = 0; testindex < 100; testindex++) {
        index query = Aux::Random::integer(n);
        if (query == n) query--;
        vector<index> lastNeighbours;
        for (double threshold = 0; threshold < R; threshold += 0.01) {
            vector<index> neighbours = quad.getElementsInHyperbolicCircle(HyperbolicSpace::polarToCartesian(angles[query], radii[query]), threshold);
            EXPECT_GE(neighbours.size(), lastNeighbours.size());
            if (neighbours.size() < lastNeighbours.size()) {
#ifndef NETWORKIT_RELEASE_LOGGING
                DEBUG("Previous Neighbours: ");
                for (index neighbour : lastNeighbours) {
                    DEBUG(neighbour);
                }
                DEBUG("Current Neighbours: ");
                for (index neighbour : neighbours) {
                    DEBUG(neighbour);
                }
#endif
            }
            lastNeighbours = neighbours;
        }
    }
}

/**
 * Insert nodes into Quadtree and successively delete all of them, check if resulting tree is empty
 */
TEST_F(QuadTreeGTest, testQuadTreeDeletion) {
    count n = 1000;
    double R = HyperbolicSpace::hyperbolicAreaToRadius(n);
    double r = HyperbolicSpace::hyperbolicRadiusToEuclidean(R);
    vector<double> angles(n);
    vector<double> radii(n);
    vector<double> indices(n);
    HyperbolicSpace::fillPoints(angles, radii, R, 1);
    double max = 0;
    for (index i = 0; i < n; i++) {
        indices[i] = i;
        radii[i] = HyperbolicSpace::hyperbolicRadiusToEuclidean(radii[i]);
        if (radii[i] == r) radii[i] = std::nextafter(radii[i], 0);

        if (radii[i] > max) {
            max = radii[i];
        }
    }

    Quadtree<index> quad(HyperbolicSpace::hyperbolicRadiusToEuclidean(R));

    for (index i = 0; i < n; i++) {
        EXPECT_GE(angles[i], 0);
        EXPECT_LT(angles[i], 2*PI);
        EXPECT_GE(radii[i], 0);
        EXPECT_LT(radii[i], R);
        TRACE("Added (", angles[i], ",", radii[i], ")");
        quad.addContent(i, angles[i], radii[i]);
    }

    while(indices.size() > 0) {
        //pick random point which is not yet deleted
        index toRemove = Aux::Random::integer(indices.size());
        if (toRemove == indices.size()) toRemove--;
        assert(toRemove < indices.size());

        //remove content at point
        EXPECT_EQ(quad.size(), indices.size());
        bool removed = quad.removeContent(indices[toRemove], angles[toRemove], radii[toRemove]);
        EXPECT_TRUE(removed);
        EXPECT_EQ(quad.size(), indices.size()-1);

        //removing twice should not work
        bool removedTwice = quad.removeContent(indices[toRemove], angles[toRemove], radii[toRemove]);
        EXPECT_FALSE(removedTwice);
        EXPECT_EQ(quad.size(), indices.size()-1);

        //mark point as removed in query list
        indices.erase(indices.begin()+toRemove);
        angles.erase(angles.begin()+toRemove);
        radii.erase(radii.begin()+toRemove);
    }

    QuadNode<index> root = getRoot(quad);
    //if root is leaf node, the coarsening worked
    EXPECT_EQ(getChildren(root).size(), 0u);
}

/**
 * Test whether the points found by a Euclidean range query on the quadtree root are exactly those whose Euclidean distance to the query point is smaller than the threshold.
 */
TEST_F(QuadTreeGTest, testEuclideanCircle) {
    // In Windows this test randomly timeouts. We temporarily address this issue
    // by fixing a random seed. However, this test should not depend on the
    // random seed.
    // TODO: debug this test on Windows.
    Aux::Random::setSeed(1, false);
    count n = 1000;
    double R = 1;
    double r = HyperbolicSpace::hyperbolicRadiusToEuclidean(R);
    vector<double> angles(n);
    vector<double> radii(n);
    HyperbolicSpace::fillPoints(angles, radii, R, 1);
    double max = 0;
    for (index i = 0; i < n; i++) {
        radii[i] = HyperbolicSpace::hyperbolicRadiusToEuclidean(radii[i]);
        if (radii[i] == r) radii[i] = std::nextafter(radii[i], 0);

        if (radii[i] > max) {
            max = radii[i];
        }
    }
    Quadtree<index> quad(max+(1-max)/4);

    for (index i = 0; i < n; i++) {
        EXPECT_GE(angles[i], 0);
        EXPECT_LT(angles[i], 2*PI);
        EXPECT_GE(radii[i], 0);
        EXPECT_LT(radii[i], R);
        TRACE("Added (", angles[i], ",", radii[i], ")");
        quad.addContent(i, angles[i], radii[i]);
    }
    vector<index> all = quad.getElements();
    EXPECT_EQ(all.size(), n);
    QuadNode<index> root = getRoot(quad);
    for (index i = 0; i < 100; i++) {
        index comparison = Aux::Random::integer(n);

        Point2DWithIndex<double> query = HyperbolicSpace::polarToCartesian(angles[comparison], radii[comparison]);
        double radius = Aux::Random::real(1);//this may overshoot the poincar disc, this is intentional. I want to know what happens
        double minR = query.length() - radius;
        double maxR = query.length() + radius;
        double minPhi, maxPhi, phi_c, r_c, spread;
        if (minR < 0) {
            maxR = std::max(abs(minR), maxR);
            minR = 0;
            minPhi = 0;
            maxPhi = 2*PI;
        } else {
            spread = asin(radius / query.length());
            HyperbolicSpace::cartesianToPolar(query, phi_c, r_c);
            minPhi = phi_c - spread;
            maxPhi = phi_c + spread;
            /**
             * what to do if they overlap the 2\pi line? Well, have to make two separate calls and collect
             */
        }

        /**
         * get Elements in Circle
         */

        vector<index> circleDenizens;

        root.getElementsInEuclideanCircle(query, radius, circleDenizens, minPhi, maxPhi, minR, maxR);
        if (minPhi < 0) {
            root.getElementsInEuclideanCircle(query, radius, circleDenizens, 2*PI+minPhi, 2*PI, minR, maxR);
        }
        if (maxPhi > 2*PI) {
            root.getElementsInEuclideanCircle(query, radius, circleDenizens, 0, maxPhi - 2*PI, minR, maxR);
        }

        //check whether bounds were correct by calling again without bounds and comparing
        vector<index> alternateDenizens;
        root.getElementsInEuclideanCircle(query, radius, alternateDenizens);

        EXPECT_EQ(circleDenizens.size(), alternateDenizens.size());


        for (index j = 0; j < n; j++) {
            Point2DWithIndex<double> comp = HyperbolicSpace::polarToCartesian(angles[j], radii[j]);
            double dist = comp.distance(query);
            if (dist < radius) {
                bool found = false;
                for (index k = 0; k < circleDenizens.size(); k++) {
                    if (circleDenizens[k] == j) {
                        found = true;
                    }
                }
                EXPECT_TRUE(found)<< "dist(" << j << "," << comparison << ") = "
                        << dist << " < " << radius;
                if (!found) {
                    DEBUG("angle: ", angles[j], ", radius: ", radii[j]);
                }
            }
        }
    }
}

/**
 * Test whether the theoretical splitting rules for different point distributions succeed in creating a balanced tree
 */
TEST_F(QuadTreeGTest, testQuadTreeBalance) {
    count n = 100000;
    double s =1;
    double alpha = 1;
    double R = s*HyperbolicSpace::hyperbolicAreaToRadius(n);
    double r = HyperbolicSpace::hyperbolicRadiusToEuclidean(R);
    vector<double> angles(n);
    vector<double> radii(n);
    HyperbolicSpace::fillPoints(angles, radii, R, alpha);
    double max = 0;
    for (index i = 0; i < n; i++) {
        radii[i] = HyperbolicSpace::hyperbolicRadiusToEuclidean(radii[i]);
        if (radii[i] == r) radii[i] = std::nextafter(radii[i], 0);
        if (radii[i] > max) {
            max = radii[i];
        }
    }

    Quadtree<index> quad(HyperbolicSpace::hyperbolicRadiusToEuclidean(R),true,alpha, 1000, 0.5);

    for (index i = 0; i < n; i++) {
        EXPECT_GE(angles[i], 0);
        EXPECT_LT(angles[i], 2*PI);
        EXPECT_GE(radii[i], 0);
        EXPECT_LT(radii[i], R);
        TRACE("Added (", angles[i], ",", radii[i], ")");
        quad.addContent(i, angles[i], radii[i]);
    }

    QuadNode<index> root = getRoot(quad);

    //visit tree
    std::stack<QuadNode<index> > toAnalyze;
    toAnalyze.push(root);
    while (!toAnalyze.empty()) {
        QuadNode<index> current = toAnalyze.top();
        toAnalyze.pop();
        if (current.height() > 1) {
            vector<QuadNode<index> > children = getChildren(current);
            EXPECT_EQ(children.size(), 4u);

            EXPECT_LE(children[0].size(), 2*children[1].size());
            EXPECT_LE(children[0].size(), 2*children[3].size());

            EXPECT_LE(children[1].size(), 2*children[0].size());
            EXPECT_LE(children[1].size(), 2*children[2].size());

            EXPECT_LE(children[2].size(), 2*children[1].size());
            EXPECT_LE(children[2].size(), 2*children[3].size());

            EXPECT_LE(children[3].size(), 2*children[2].size());
            EXPECT_LE(children[3].size(), 2*children[0].size());

            EXPECT_LE(children[0].height(), children[1].height()+1);
            EXPECT_LE(children[0].height(), children[3].height()+1);

            EXPECT_LE(children[1].height(), children[0].height()+1);
            EXPECT_LE(children[1].height(), children[2].height()+1);

            EXPECT_LE(children[2].height(), children[1].height()+1);
            EXPECT_LE(children[2].height(), children[3].height()+1);

            EXPECT_LE(children[3].height(), children[2].height()+1);
            EXPECT_LE(children[3].height(), children[0].height()+1);
            for (auto child : children) {
                toAnalyze.push(child);
            }
        }
    }
}

TEST_F(QuadTreeGTest, testSequentialQuadTreeConstruction) {
    count n = 100000;
    count capacity = 1000;
    double alpha = 1;
    double R = HyperbolicSpace::hyperbolicAreaToRadius(n);
    double r = HyperbolicSpace::hyperbolicRadiusToEuclidean(R);
    vector<double> angles(n);
    vector<double> radii(n);
    HyperbolicSpace::fillPoints(angles, radii, R, alpha);

    Quadtree<index> quad(HyperbolicSpace::hyperbolicRadiusToEuclidean(R),false,alpha,capacity);

    for (index i = 0; i < n; i++) {
        radii[i] = HyperbolicSpace::hyperbolicRadiusToEuclidean(radii[i]);
        if (radii[i] == r) radii[i] = std::nextafter(radii[i], 0);
        quad.addContent(i, angles[i], radii[i]);
    }
    EXPECT_EQ(quad.size(), n);

    quad.trim();
    vector<double> anglecopy;
    vector<double> radiicopy;
    quad.extractCoordinates(anglecopy, radiicopy);
    EXPECT_EQ(anglecopy.size(), n);
    EXPECT_EQ(radiicopy.size(), n);
    //EXPECT_TRUE(std::is_permutation(angles.begin(), angles.end(), anglecopy.begin()));//too slow!
    //EXPECT_TRUE(std::is_permutation(radii.begin(), radii.end(), radiicopy.begin()));
}

TEST_F(QuadTreeGTest, testProbabilisticQuery) {
    Aux::Random::setSeed(0, false);
    count n = 5000;
    count m = n*3;
    count capacity = 20;
    double R = 2*log(8*n / (PI*(m/n)*2));
    double r = HyperbolicSpace::hyperbolicRadiusToEuclidean(R);
    double alpha = 1;

    vector<double> angles(n);
    vector<double> radii(n);

    HyperbolicSpace::fillPoints(angles, radii, R, alpha);

    Quadtree<index> quad(r,false,alpha,capacity);

    for (index i = 0; i < n; i++) {
        radii[i] = HyperbolicSpace::hyperbolicRadiusToEuclidean(radii[i]);
        if (radii[i] == r) radii[i] = std::nextafter(radii[i], 0);
        EXPECT_EQ(i, quad.size());
        quad.addContent(i, angles[i], radii[i]);
    }
    EXPECT_EQ(n, quad.size());

    quad.trim();

    for (index i = 0; i < 200; i++) {
        index query = Aux::Random::integer(n-1);
        double acc = Aux::Random::probability() ;
        auto edgeProb = [acc](double) -> double {return acc;};
        vector<index> near;
        quad.getElementsProbabilistically(HyperbolicSpace::polarToCartesian(angles[query], radii[query]), edgeProb, near);
        EXPECT_NEAR(near.size(), acc*n, std::max(acc*n*0.25, 10.0));
    }

    //TODO: some test about appropriate subtrees and leaves

    auto edgeProb = [](double) -> double {return 1;};
    vector<index> near;
    quad.getElementsProbabilistically(HyperbolicSpace::polarToCartesian(angles[0], radii[0]), edgeProb, near);
    EXPECT_EQ(n, near.size());

    auto edgeProb2 = [](double) -> double {return 0;};
    near.clear();
    quad.getElementsProbabilistically(HyperbolicSpace::polarToCartesian(angles[0], radii[0]), edgeProb2, near);
    EXPECT_EQ(0u, near.size());
}

TEST_F(QuadTreeGTest, testCartesianEuclidQuery) {
    count n = 10000;
    //count m = n*3;
    //count capacity = 20;

    assert(n > 0);

    vector<Point<double> > positions(n);
    vector<index> content(n);

    for (index i = 0; i < n; i++) {
        Point<double> pos{Aux::Random::probability(), Aux::Random::probability()};
        positions[i] = pos;
        content[i] = i;
    }

    QuadtreeCartesianEuclid<index> quad(positions, content, true);

    EXPECT_EQ(n, quad.size());
    quad.recount();
    EXPECT_EQ(n, quad.size());

    quad.trim();

    for (index i = 0; i < 200; i++) {
        index query = Aux::Random::integer(n-1);
        double threshold = Aux::Random::real(0, 1);
        auto edgeProb = [threshold](double distance) -> double {return distance <= threshold ? 1 : 0;};
        vector<index> near;
        quad.getElementsProbabilistically(positions[query], edgeProb, near);
        vector<index> circleDenizens;
        quad.getElementsInEuclideanCircle(positions[query], threshold, circleDenizens);
        EXPECT_EQ(near.size(), circleDenizens.size());
    }

    //TODO: some test about appropriate subtrees and leaves

    auto edgeProb = [](double) -> double {return 1;};
    vector<index> near;
    quad.getElementsProbabilistically(positions[0], edgeProb, near);
    EXPECT_EQ(n, near.size());

    auto edgeProb2 = [](double) -> double {return 0;};
    near.clear();
    quad.getElementsProbabilistically(positions[0], edgeProb2, near);
    EXPECT_EQ(0u, near.size());
}



TEST_F(QuadTreeGTest, testLeftSuppression) {
    /**
     * define parameters
     */
    count n = 10000;
    double k = 10;
    count m = n*k/2;
    double targetR = HyperbolicSpace::getTargetRadius(n, m);
    double R = HyperbolicSpace::hyperbolicAreaToRadius(n);
    double r = HyperbolicSpace::hyperbolicRadiusToEuclidean(R);
    double alpha = 1;

    //allocate data structures
    vector<double> angles(n), radii(n);

    /**
     * generate values and construct quadtree
     */
    HyperbolicSpace::fillPoints(angles, radii, targetR, alpha);
    Quadtree<index> quad(HyperbolicSpace::hyperbolicRadiusToEuclidean(targetR));

    for (index i = 0; i < n; i++) {
        radii[i] = HyperbolicSpace::hyperbolicRadiusToEuclidean(radii[i]);
        if (radii[i] == r) radii[i] = std::nextafter(radii[i], 0);
        quad.addContent(i, angles[i], radii[i]);
    }
    EXPECT_EQ(quad.size(), n);

    //now test
    for (index i = 0; i < n; i++) {
        vector<index> allNeighbours;
        vector<index> rightNeighbours;
        quad.getElementsInHyperbolicCircle(HyperbolicSpace::polarToCartesian(angles[i], radii[i]), targetR, true, rightNeighbours);
        quad.getElementsInHyperbolicCircle(HyperbolicSpace::polarToCartesian(angles[i], radii[i]), targetR, false, allNeighbours);
        EXPECT_LE(rightNeighbours.size(), allNeighbours.size());
        index aIndex, bIndex;
        for (aIndex = 0, bIndex = 0; aIndex < rightNeighbours.size() && bIndex < allNeighbours.size(); aIndex++, bIndex++)  {
            //EXPECT_GE(angles[rightNeighbours[aIndex]], angles[i]);//all elements returned by partial query are right
            while(rightNeighbours[aIndex] != allNeighbours[bIndex]) {//iterate over suppressed elements until next match
                EXPECT_LT(angles[allNeighbours[bIndex]], angles[i]);//all elements suppressed are left
                bIndex++;
            }
        }
        EXPECT_EQ(aIndex, rightNeighbours.size());//all elements in partial query are also contained in complete query
    }
}

TEST_F(QuadTreeGTest, debugTreeExport) {
    count n = 200;
    count capacity = 40;
    double k = 10;
    count m = n*k/2;
    double targetR = HyperbolicSpace::getTargetRadius(n, m);

    double R = HyperbolicSpace::hyperbolicAreaToRadius(n);
    double r = HyperbolicSpace::hyperbolicRadiusToEuclidean(R);
    double alpha = 0.7;

    //allocate data structures
    vector<double> angles(n), radii(n);

    /**
     * generate values and construct quadtree
     */
    HyperbolicSpace::fillPoints(angles, radii, targetR, alpha);
    Quadtree<index> quad(HyperbolicSpace::hyperbolicRadiusToEuclidean(targetR), true, alpha, capacity);

    for (index i = 0; i < n; i++) {
        radii[i] = HyperbolicSpace::hyperbolicRadiusToEuclidean(radii[i]);
        if (radii[i] == r) radii[i] = std::nextafter(radii[i], 0);
        quad.addContent(i, angles[i], radii[i]);
    }

    EXPECT_EQ(quad.size(), n);

    quad.indexSubtree(0);

    count treeheight = quad.height();
    DEBUG("Quadtree height: ", treeheight);
    auto deg = [](double rad) -> double {return 180*rad/PI;};
    tlx::unused(treeheight, deg);

    index query = Aux::Random::integer(n-1);
    radii[query] = HyperbolicSpace::hyperbolicRadiusToEuclidean(HyperbolicSpace::EuclideanRadiusToHyperbolic(radii[query])*0.75);
    DEBUG("Query:", angles[query], ", ", radii[query]);
    double T = 0.5;
    double beta = 1/T;

    auto edgeProb = [beta, targetR](double distance) -> double {return 1 / (exp(beta*(distance-targetR)/2)+1);};

    std::stack<std::tuple<QuadNode<index>, count, double, double, index > > quadnodestack;
    QuadNode<index> root = getRoot(quad);
    quadnodestack.push(std::make_tuple(root, quad.height(), 0, 0, root.getID()));

    while(!quadnodestack.empty()) {
        auto currentTuple = quadnodestack.top();
        quadnodestack.pop();
        QuadNode<index> current;
        count remainingHeight;
        double xoffset, yoffset;
        index parentID;
        std::tie(current, remainingHeight, xoffset, yoffset, parentID) = currentTuple;

        DEBUG("Quadtree Cell ", current.getID());
        double probUB = edgeProb(current.hyperbolicDistances(angles[query], radii[query]).first);
        DEBUG("\\drawQuadNode{",xoffset,"}{", yoffset, "}{", current.getID(), "}{", parentID, "}{", probUB, "}{", current.size(), "}");
        DEBUG("\\drawCell{",deg(current.getLeftAngle()), "}{", deg(current.getRightAngle()),"}{", HyperbolicSpace::EuclideanRadiusToHyperbolic(current.getMinR()), "}{", HyperbolicSpace::EuclideanRadiusToHyperbolic(current.getMaxR()), "}");
        DEBUG("Height: ", current.height());
        auto distances = current.hyperbolicDistances(angles[query], radii[query]);
        DEBUG("Mindistance to query:", distances.first);
        DEBUG("ProbUB:", edgeProb(distances.first), " ProbLB:", edgeProb(distances.second));
        tlx::unused(probUB, distances);

        if (current.height() == 1) {
            index i = 0;
            for (index elem : current.getElements()) {
                i++;
                double p = edgeProb(HyperbolicSpace::poincareMetric(angles[elem], radii[elem], angles[query], radii[query]));
                DEBUG("\\drawPoint{", deg(angles[elem]), "}{", HyperbolicSpace::EuclideanRadiusToHyperbolic(radii[elem]), "}{", p, "}{", current.getID(), "}{", i, "}");
                DEBUG("Leaf contains: ", angles[elem], ", ", radii[elem], " p: ", p);
                tlx::unused(p);
            }
        }

        double stepsize = pow(4, remainingHeight-1);
        double newXOffset = xoffset-1.5*stepsize;
        double newYOffset = yoffset + 1;
        for (QuadNode<index> child : current.children) {
            quadnodestack.push(std::make_tuple(child, remainingHeight-1, newXOffset, newYOffset, current.getID()));
            newXOffset += stepsize;
        }
    }


    DEBUG("\\drawQuery{", deg(angles[query]), "}{", HyperbolicSpace::EuclideanRadiusToHyperbolic(radii[query]),"}");

}

TEST_F(QuadTreeGTest, testPolarEuclidQuery) {
    Aux::Random::setSeed(42, true);
    /**
     * setup of data structures and constants
     */
    double maxR = 2;
    count n = 10000;
    vector<double> angles(n);
    vector<double> radii(n);
    vector<index> content(n);

    double minPhi = 0;
    double maxPhi = 2*PI;
    double minR = 0;

    /**
     * get random number generators
     */

    std::uniform_real_distribution<double> phidist{minPhi, maxPhi};
    std::uniform_real_distribution<double> rdist{minR, maxR};

    /**
     * fill vectors
     */
    for (index i = 0; i < n; i++) {
        angles[i] = phidist(Aux::Random::getURNG());
        radii[i] = rdist(Aux::Random::getURNG());
        content[i] = i;
    }

    const bool splitTheoretical = true;
    QuadtreePolarEuclid<index> tree(angles, radii, content, splitTheoretical);
    EXPECT_EQ(n, tree.size());

    tree.trim();

    for (index i = 0; i < 200; i++) {
        index query = Aux::Random::integer(n-1);
        double acc = Aux::Random::probability() ;
        auto edgeProb = [acc](double) -> double {return acc;};
        vector<index> near;
        tree.getElementsProbabilistically(HyperbolicSpace::polarToCartesian(angles[query], radii[query]), edgeProb, near);
        EXPECT_NEAR(near.size(), acc*n, std::max(acc*n*0.25, 10.0));
    }

    //TODO: some test about appropriate subtrees and leaves

    auto edgeProb = [](double) -> double {return 1;};
    vector<index> near;
    tree.getElementsProbabilistically(HyperbolicSpace::polarToCartesian(angles[0], radii[0]), edgeProb, near);
    EXPECT_EQ(n, near.size());

    auto edgeProb2 = [](double) -> double {return 0;};
    near.clear();
    tree.getElementsProbabilistically(HyperbolicSpace::polarToCartesian(angles[0], radii[0]), edgeProb2, near);
    EXPECT_EQ(0u, near.size());
}

TEST_F(QuadTreeGTest, testQuadTreePolarEuclidInsertion) {
    /**
     * setup of data structures and constants
     */
    double maxR = 2;
    count n = 1000;
    vector<double> angles(n);
    vector<double> radii(n);
    vector<index> content(n);

    double minPhi = 0;
    double maxPhi = 2*PI;
    double minR = 0;

    /**
     * get random number generators
     */

    std::uniform_real_distribution<double> phidist{minPhi, maxPhi};
    std::uniform_real_distribution<double> rdist{minR, maxR};

    /**
     * fill vectors
     */
    for (index i = 0; i < n; i++) {
        angles[i] = phidist(Aux::Random::getURNG());
        radii[i] = rdist(Aux::Random::getURNG());
        content[i] = i;
    }

    QuadtreePolarEuclid<index> tree(angles, radii, content);
    EXPECT_EQ(n, tree.size());

    /**
     * elements are returned
     */
    vector<index> returned = tree.getElements();
    EXPECT_EQ(n, returned.size());
    sort(returned.begin(), returned.end());
    for (index i = 0; i < returned.size(); i++) {
        EXPECT_EQ(i, returned[i]);
    }
}

TEST_F(QuadTreeGTest, testQuadNodeHyperbolicDistances) {
    double minR, maxR, minPhi, maxPhi;
    std::tie(minPhi, minR, maxPhi, maxR) = std::make_tuple(6.1850105367549055, 0.0, 6.2831853071795862, 0.99998090224083802);
    ASSERT_LE(minPhi, maxPhi);
    ASSERT_LE(minR, maxR);

    Point2DWithIndex<double> query(0.99999879353779952, 0.0011511619164281455);
    double phi_q, r_q;
    HyperbolicSpace::cartesianToPolar(query, phi_q, r_q);

    QuadNode<index> node(minPhi, minR, maxPhi, maxR, 1000, 0);
    count steps = 100;
    Point2DWithIndex<double> posAtMin = HyperbolicSpace::polarToCartesian(minPhi, minR);
    double minDistance = HyperbolicSpace::poincareMetric(query, posAtMin);
    double minR_Hyper = HyperbolicSpace::EuclideanRadiusToHyperbolic(minR);

    double phiStep = (maxPhi-minPhi)/steps;
    double rStep = (HyperbolicSpace::EuclideanRadiusToHyperbolic(maxR) - minR_Hyper) / steps;
    for (index i = 0; i <= steps; i++) {
        double phi = minPhi + i*phiStep;
        for (index j = 0; j <= steps; j++) {
            double r =  HyperbolicSpace::hyperbolicRadiusToEuclidean(minR_Hyper + j*rStep);
            if (i < steps && j < steps) {
                EXPECT_TRUE(node.responsible(phi, r));
            } else {
                EXPECT_FALSE(node.responsible(phi, r));
            }
            Point2DWithIndex<double> pos = HyperbolicSpace::polarToCartesian(phi, r);
            double distance =  HyperbolicSpace::poincareMetric(query, pos);
            if (distance < minDistance) {
                minDistance = distance;
                posAtMin = pos;
            }
        }
    }

    double phiAtMin, rAtMin;
    HyperbolicSpace::cartesianToPolar(posAtMin, phiAtMin, rAtMin);
    DEBUG("Point in Cell at minimal distance at (", phiAtMin, ", ", rAtMin, ") on the poincare disk, which is ", HyperbolicSpace::EuclideanRadiusToHyperbolic(rAtMin), " in native radius.");

    Point2DWithIndex<double> p(0.99983582188685627, -0.0020919222764388666);
    double phi, r;
    HyperbolicSpace::cartesianToPolar(p, phi, r);

    EXPECT_TRUE(node.responsible(phi, r));
    double distanceQueryToCell = node.hyperbolicDistances(phi_q, r_q).first;
    double distanceQueryToPoint = HyperbolicSpace::poincareMetric(query, p);

    //node.addContent(1, phi, r);
    EXPECT_LE(distanceQueryToCell, distanceQueryToPoint);


}

TEST_F(QuadTreeGTest, testQuadNodeCartesianDistances) {
    Point<double> lower(0.24997519780061023, 0.7499644402803205);
    Point<double> upper(0.49995039560122045, 0.99995258704042733);

    ASSERT_LE(lower[0], upper[0]);
    ASSERT_LE(lower[1], upper[1]);
    ASSERT_EQ(2u, lower.getDimensions());
    ASSERT_EQ(2u, upper.getDimensions());

    Point<double> query(0.81847946542324035, 0.91885035291473593);

    QuadNodeCartesianEuclid<index> node(lower, upper, 1000);
    //count steps = 100;
    Point<double> posAtMin = lower;
    double minDistance = posAtMin.distance(query);

    DEBUG("Point in Cell at minimal distance at (", posAtMin[0], ", ", posAtMin[0], "), distance is ", minDistance);

    Point<double> p(0.49969783875749996, 0.87199796797360407);

    EXPECT_TRUE(node.responsible(p));
    double distanceQueryToCell = node.EuclideanDistances(query).first;
    double distanceQueryToPoint = query.distance(p);

    EXPECT_LE(distanceQueryToCell, distanceQueryToPoint);
    EXPECT_LE(distanceQueryToCell, minDistance);


}

} /* namespace NetworKit */
