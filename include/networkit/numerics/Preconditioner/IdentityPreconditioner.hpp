// no-networkit-format
/*
 * IdentityPreconditioner.h
 *
 *  Created on: Apr 23, 2016
 *      Author: Michael Wegner (michael.wegner@student.kit.edu)
 */

#ifndef NETWORKIT_NUMERICS_PRECONDITIONER_IDENTITY_PRECONDITIONER_HPP_
#define NETWORKIT_NUMERICS_PRECONDITIONER_IDENTITY_PRECONDITIONER_HPP_

namespace NetworKit {

/**
 * @ingroup numerics
 * Simple preconditioner that returns the given vector unchanged.
 */
class IdentityPreconditioner {
public:
    /** Default constructor */
    IdentityPreconditioner() = default;
    /**
     * Constructs an identity preconditioner for the matrix @a A.
     * @param A
     */
    IdentityPreconditioner(const CSRMatrix &)  {}
    virtual ~IdentityPreconditioner() = default;

    /**
     * Returns the preconditioned right-hand-side \f$P(b) = b\f$.
     */
    Vector rhs(const Vector& b) const {
        return b;
    }
};

} /* namespace NetworKit */

#endif // NETWORKIT_NUMERICS_PRECONDITIONER_IDENTITY_PRECONDITIONER_HPP_
