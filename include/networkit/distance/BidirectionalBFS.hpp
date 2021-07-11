// no-networkit-format
/*
 * BidirectionalBFS.hpp
 *
 *  Created on: 14.06.2019
 *      Author: Eugenio Angriman <angrimae@hu-berlin.de>
 */

#ifndef NETWORKIT_DISTANCE_BIDIRECTIONAL_BFS_HPP_
#define NETWORKIT_DISTANCE_BIDIRECTIONAL_BFS_HPP_

#include <queue>

#include <networkit/distance/STSP.hpp>

namespace NetworKit {

/*
 * @ingroup distance
 * Implements a bidirectional breadth-first search on a graph from two given source and target
 * nodes. Explores the graph from both the source and target nodes until the two explorations meet.
 */
class BidirectionalBFS final : public STSP {

public:
    // Inherit the constructors of STSP.
    using STSP::STSP;

    /*
     * Perform a bidirectional BFS from the given source and target nodes.
     */
    void run() override;

    /*
     * Returns the distance (i.e., number of hops) from the source to the target node.
     *
     * @return count Number of hops from the source to the target node.
     */
    count getHops() {
        assureFinished();
        return stDist;
    }

    edgeweight getDistance() const override {
        assureFinished();
        return static_cast<edgeweight>(stDist);
    }

private:
    count stDist;
    std::vector<uint8_t> visited;
    uint8_t ts = 0;
    static constexpr uint8_t ballMask = uint8_t(1) << 7;
    std::queue<node> sQueue, sQueueNext, tQueue, tQueueNext;
};
} // namespace NetworKit

#endif // NETWORKIT_DISTANCE_BIDIRECTIONAL_BFS_HPP_
