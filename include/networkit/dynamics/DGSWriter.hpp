// no-networkit-format
/*
 * DGSWriter.hpp
 *
 *  Created on: 14.01.2014
 *      Author: cls
 */

#ifndef NETWORKIT_DYNAMICS_DGS_WRITER_HPP_
#define NETWORKIT_DYNAMICS_DGS_WRITER_HPP_

#include <string>
#include <vector>

#include <networkit/dynamics/GraphEvent.hpp>

namespace NetworKit {

/**
 * @ingroup dynamics
 */
class DGSWriter final {
public:
    DGSWriter() = default;

    void write(std::vector<GraphEvent>& stream, const std::string& path);
};

} /* namespace NetworKit */

#endif // NETWORKIT_DYNAMICS_DGS_WRITER_HPP_
