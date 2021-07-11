// no-networkit-format
/*
 * DynSSSP.cpp
 *
 *  Created on: 17.07.2014
 *      Author: cls, ebergamini
 */


#include <networkit/distance/DynSSSP.hpp>

namespace NetworKit {

    DynSSSP::DynSSSP(const Graph& G, node source, bool storePredecessors, node target) : SSSP(G, source, true, false, target),
    storePreds(storePredecessors) {

    }

} /* namespace NetworKit */
