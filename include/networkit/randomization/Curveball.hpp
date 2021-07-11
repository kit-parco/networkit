/*
 * Curveball.hpp
 *
 *  Created on: 26.05.2018
 *      Author:  Hung Tran <htran@ae.cs.uni-frankfurt.de>, Manuel Penschuck <networkit@manuel.jetzt>
 */

#ifndef NETWORKIT_RANDOMIZATION_CURVEBALL_HPP_
#define NETWORKIT_RANDOMIZATION_CURVEBALL_HPP_

#include <memory>
#include <utility>

#include <networkit/Globals.hpp>
#include <networkit/base/Algorithm.hpp>
#include <networkit/graph/Graph.hpp>

#include <tlx/define/deprecated.hpp>

namespace NetworKit {

// pImpl
namespace CurveballDetails {
class CurveballIM;
}

class Curveball final : public Algorithm {
public:
    explicit Curveball(const Graph &G);

    ~Curveball() override;

    void run() final {
        throw std::runtime_error("run() is not supported by this algorithm; use run(trades)");
    };

    /* The algorithm is not parallel. */
    void run(const std::vector<std::pair<node, node>> &trades);

    Graph getGraph(bool parallel = false);

    std::string TLX_DEPRECATED(toString() const final);

    bool TLX_DEPRECATED(isParallel() const final) { return false; }

    count getNumberOfAffectedEdges() const;

private:
    std::unique_ptr<CurveballDetails::CurveballIM> impl;
};

} // namespace NetworKit

#endif // NETWORKIT_RANDOMIZATION_CURVEBALL_HPP_
