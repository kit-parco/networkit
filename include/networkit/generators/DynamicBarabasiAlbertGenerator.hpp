// no-networkit-format
/*
 * DynamicBarabasiAlbertGenerator.hpp
 *
 *  Created on: 09.04.2013
 *      Author: cls
 */

#ifndef NETWORKIT_GENERATORS_DYNAMIC_BARABASI_ALBERT_GENERATOR_HPP_
#define NETWORKIT_GENERATORS_DYNAMIC_BARABASI_ALBERT_GENERATOR_HPP_

#include <set>

#include <networkit/generators/DynamicGraphSource.hpp>

namespace NetworKit {


// FIXME: for k=2, degree 2 nodes should be most frequent but degree 4 nodes are
/**
 * @ingroup generators
 */
class DynamicBarabasiAlbertGenerator final : public DynamicGraphSource {

    count k;             //!< parameter of the BA model: number of edges per new node
    count degSum;        //!< degree sum of current graph

public:

    DynamicBarabasiAlbertGenerator(count k = 2);

    void initializeGraph() override;

    void generate() override;

};

} /* namespace NetworKit */
#endif // NETWORKIT_GENERATORS_DYNAMIC_BARABASI_ALBERT_GENERATOR_HPP_
