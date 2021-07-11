// no-networkit-format
/*
 * DynamicPubWebGenerator.cpp
 *
 *  Created on: 15.01.2014
 *      Author: Henning
 */

#include <map>

#include <networkit/generators/DynamicPubWebGenerator.hpp>

namespace NetworKit {

DynamicPubWebGenerator::DynamicPubWebGenerator(count numNodes, count numberOfDenseAreas,
                                               coordinate neighborhoodRadius, count maxNumberOfNeighbors,
                                               bool writeInitialGraphToStream)
    : initGen(numNodes, numberOfDenseAreas, neighborhoodRadius, maxNumberOfNeighbors),
      writeInitialGraphToStream(writeInitialGraphToStream) {
    G = initGen.generate();
    coordinates = initGen.moveCoordinates();
}

std::vector<GraphEvent> DynamicPubWebGenerator::generate(count nSteps) {
    const auto numToDel = static_cast<count>(G.numberOfNodes() * 0.05);
    const auto numToIns = static_cast<count>(G.numberOfNodes() * 0.05);
    std::vector<GraphEvent> eventStream;
    std::vector<std::pair<node, node>> edgesToDelete;
    newCoordinates.clear();
    newCoordinates.reserve(numToIns);

    if (writeInitialGraphToStream) {
        // write initial graph to stream
        G.forNodes([&](node v) { eventStream.emplace_back(GraphEvent::NODE_ADDITION, v); });

        G.forEdges([&](node u, node v, edgeweight ew) {
            eventStream.emplace_back(GraphEvent::EDGE_ADDITION, u, v, ew);
        });

        eventStream.emplace_back(GraphEvent::TIME_STEP);

        writeInitialGraphToStream = false;
    }

    for (index step = 0; step < nSteps; ++step) {

        // delete nodes
        for (index i = 0; i < numToDel; ++i) {
            // draw a certain (random) number of vertices to be deleted
            node nodeToDel = none;

            do {
                nodeToDel = Aux::Random::integer(G.upperNodeIdBound() - 1);
            } while (!G.hasNode(nodeToDel));

            // mark incident edges first for deletion, delete them outside of neighborhood iterator
            G.forNeighborsOf(nodeToDel,
                             [&](node neigh) { edgesToDelete.emplace_back(nodeToDel, neigh); });
            for (auto elem : edgesToDelete) {
                node neigh = elem.second;
                G.removeEdge(nodeToDel, neigh);
                GraphEvent event(GraphEvent::EDGE_REMOVAL, nodeToDel, neigh);
                eventStream.push_back(event);
            }
            edgesToDelete.clear();

            // eventually delete vertex
            G.removeNode(nodeToDel);
            GraphEvent event(GraphEvent::NODE_REMOVAL, nodeToDel);
            eventStream.push_back(event);
        }

        // insert nodes
        coordinates.reserve(coordinates.size() + numToIns);
        for (index i = 0; i < numToIns; ++i) {
            // draw a cluster where the vertex should be inserted, +1 to account for the noise
            auto drawCoordinate = [&]() -> Point2D {
                count clusterToIns = Aux::Random::integer(initGen.numDenseAreas + 1);

                if (clusterToIns < initGen.numDenseAreas) {
                    // real cluster, FIXME: DRY!
                    // compute random angle between [0, 2pi) and distance between [0, width/2]
                    coordinate angle = Aux::Random::real() * 2.0 * PI;
                    coordinate dist = Aux::Random::real() * initGen.denseAreaXYR[clusterToIns].rad;

                    // compute coordinates and adjust them
                    return initGen.intoUnitSquare(
                        {initGen.denseAreaXYR[clusterToIns].x + std::cos(angle) * dist,
                         initGen.denseAreaXYR[clusterToIns].y + std::sin(angle) * dist});
                } else {
                    // noise -> random coordinate
                    return {Aux::Random::probability(), Aux::Random::probability()};
                }
            };

            // create vertex with these coordinates
            node newNode = G.addNode();
            const auto coord = drawCoordinate();
            coordinates.emplace_back(coord);
            newCoordinates.emplace_back(newNode, coord);
            GraphEvent event(GraphEvent::NODE_ADDITION, newNode);
            eventStream.push_back(event);
        }

        // determine events by computing new graph structure
        using edge = std::pair<node, node>;
        coordinate sqrNeighRad = initGen.neighRad * initGen.neighRad;
        std::map<edge, count> eligibleEdges;

        // find for each node the rad-neighborhood
        // FIXME: get rid of quadratic running time!
        G.forNodes([&](node u) {
            std::priority_queue<std::pair<coordinate, edge>> pq;

            // fill PQ with neighbors in range
            G.forNodes([&](node v) {
                coordinate sqrDist = initGen.squaredDistanceInUnitTorus(coordinates[u], coordinates[v]);

                if (sqrDist <= sqrNeighRad) {
                    edge e{std::min(u, v), std::max(u, v)};
                    pq.emplace(-sqrDist, e);
                }
            });

            // mark up to maxNeigh nearest neighbors as eligible
            count end = std::min(initGen.maxNeigh, (count)pq.size());
            for (index i = 0; i < end; ++i) {
                eligibleEdges[pq.top().second]++;
                pq.pop();
            }
        });

        // check if edges have to be deleted (was there, but not eligible twice any more)
        G.forEdges([&](node u, node v) {
            edge e = std::make_pair(std::min(u, v), std::max(u, v));
            if (eligibleEdges[e] < 2) {
                edgesToDelete.emplace_back(u, v);
            } else {
                assert(G.hasEdge(u, v));
                edgeweight ew =
                    PubWebGenerator::BASE_WEIGHT
                    / initGen.squaredDistanceInUnitTorus(coordinates[u], coordinates[v]);
                G.setWeight(u, v, ew);
                GraphEvent event(GraphEvent::EDGE_WEIGHT_UPDATE, u, v, ew);
                eventStream.push_back(event);
            }
            eligibleEdges.erase(e);
        });

        // now really delete the edges after leaving the G.forEdges iterator
        for (auto elem : edgesToDelete) {
            node &u = elem.first;
            node &v = elem.second;
            G.removeEdge(u, v);
            GraphEvent event(GraphEvent::EDGE_REMOVAL, u, v);
            eventStream.push_back(event);
        }
        edgesToDelete.clear();

        // check if edges have to be inserted
        for (const auto &edgePair : eligibleEdges) {
            if (edgePair.second >= 2) {
                node u, v;
                std::tie(u, v) = edgePair.first;
                assert(!G.hasEdge(u, v) && G.hasNode(u) && G.hasNode(v));

                edgeweight ew =
                    PubWebGenerator::BASE_WEIGHT
                    / initGen.squaredDistanceInUnitTorus(coordinates[u], coordinates[v]);
                G.addEdge(u, v, ew);
                assert(G.hasEdge(u, v));
                GraphEvent event(GraphEvent::EDGE_ADDITION, u, v, ew);
                eventStream.push_back(event);
            }
        }

        eventStream.emplace_back(GraphEvent::TIME_STEP);
    }

    return eventStream;
}

} /* namespace NetworKit */
