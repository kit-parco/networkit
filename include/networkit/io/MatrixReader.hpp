/*
 * MatrixReader.hpp
 *
 *  Created on: 17.01.2013
 *      Author: Christian Staudt
 */

#ifndef NETWORKIT_IO_MATRIX_READER_HPP_
#define NETWORKIT_IO_MATRIX_READER_HPP_

#include <networkit/algebraic/CSRMatrix.hpp>

namespace NetworKit {

/**
 * @ingroup io
 * Abstract base class for matrix readers.
 */
class MatrixReader {
public:
    virtual ~MatrixReader() = default;

    /**
     * Reads the matrix in @a path.
     */
    virtual CSRMatrix read(const std::string &path) = 0;

    /** only to be used by Cython - this eliminates an unnecessary copy */
    CSRMatrix *_read(const std::string &path) { return new CSRMatrix{read(path)}; };
};

} // namespace NetworKit
#endif // NETWORKIT_IO_MATRIX_READER_HPP_
