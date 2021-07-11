// no-networkit-format
/*
 * StaticDegreeSequenceGenerator.hpp
 *
 *  Created on: 24.02.2014
 *      Author: Henning
 */

#ifndef NETWORKIT_GENERATORS_STATIC_DEGREE_SEQUENCE_GENERATOR_HPP_
#define NETWORKIT_GENERATORS_STATIC_DEGREE_SEQUENCE_GENERATOR_HPP_

#include <networkit/generators/StaticGraphGenerator.hpp>

namespace NetworKit {

// TODO: Clean this up.
const short NO = 0;
const short YES = 1;
const short UNKNOWN = 2;

/**
 * @ingroup generators
 */
class StaticDegreeSequenceGenerator: public StaticGraphGenerator {
protected:
    std::vector<count> seq;
    short realizable;


public:
    StaticDegreeSequenceGenerator(const std::vector<count>& sequence);

    /**
     * Erdoes-Gallai test if degree sequence seq is realizable.
     */
    virtual bool isRealizable();

    virtual bool getRealizable() const;

    Graph generate() override = 0;
};

} /* namespace NetworKit */
#endif // NETWORKIT_GENERATORS_STATIC_DEGREE_SEQUENCE_GENERATOR_HPP_
