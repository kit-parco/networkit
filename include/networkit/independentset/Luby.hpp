// no-networkit-format
/*
 * Luby.hpp
 *
 *  Created on: 27.02.2013
 *      Author: Christian Staudt
 */

#ifndef NETWORKIT_INDEPENDENTSET_LUBY_HPP_
#define NETWORKIT_INDEPENDENTSET_LUBY_HPP_

#include <networkit/independentset/IndependentSetFinder.hpp>


namespace NetworKit {

/**
 * @ingroup independentset
 *
 * Luby's parallel independent set algorithm.
 */
class Luby final : public IndependentSetFinder {

public:

    // FIXME: check correctness of implementation
    std::vector<bool> run(const Graph& G) override;

    std::string TLX_DEPRECATED(toString() const override);
};

} /* namespace NetworKit */
#endif // NETWORKIT_INDEPENDENTSET_LUBY_HPP_
