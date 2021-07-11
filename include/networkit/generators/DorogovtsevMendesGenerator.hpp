// no-networkit-format
/*
* DorogovtsevMendesGenerator.hpp
*
*  Created on: 27.05.2014
*      Author: Christian Staudt
*/

#ifndef NETWORKIT_GENERATORS_DOROGOVTSEV_MENDES_GENERATOR_HPP_
#define NETWORKIT_GENERATORS_DOROGOVTSEV_MENDES_GENERATOR_HPP_

#include <networkit/generators/StaticGraphGenerator.hpp>

namespace NetworKit {

/**
 * @ingroup generators
 */
class DorogovtsevMendesGenerator final : public StaticGraphGenerator {

public:
    /**
    * @param nNodes number of nodes in target graph
    */
    DorogovtsevMendesGenerator(count nNodes);

    Graph generate() override;

private:
        count nNodes;

};

} /* namespace NetworKit */
#endif // NETWORKIT_GENERATORS_DOROGOVTSEV_MENDES_GENERATOR_HPP_
