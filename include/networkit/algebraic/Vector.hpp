/*
 * Vector.hpp
 *
 *  Created on: 12.03.2014
 *      Author: Michael Wegner (michael.wegner@student.kit.edu)
 */

#ifndef NETWORKIT_ALGEBRAIC_VECTOR_HPP_
#define NETWORKIT_ALGEBRAIC_VECTOR_HPP_

#include <cassert>
#include <cmath>
#include <stdexcept>
#include <vector>

#include <networkit/Globals.hpp>
#include <networkit/algebraic/AlgebraicGlobals.hpp>

namespace NetworKit {

// forward declaration of DynamicMatrix class
class DynamicMatrix;

/**
 * @ingroup algebraic
 * The Vector class represents a basic vector with double coefficients.
 */
class Vector final {
private:
    std::vector<double> values;
    bool transposed;

public:
    /** Default constructor */
    Vector();

    /**
     * Constructs the Vector with @a dimension elements with value @a initialValue.
     * @param dimension The dimension of this vector.
     * @param initialValue All coefficients will be initialized to @a initialValue.
     * @param transpose Indicates whether this vector is transposed (row vector) or not (column
     * vector).
     */
    Vector(count dimension, double initialValue = 0, bool transpose = false);

    /**
     * Constructs the Vector with the contents of @a values.
     * @param values The values of this Vector.
     * @param transpose Indicates whether this vector is transposed (row vector) or not (column
     * vector).
     */
    Vector(const std::vector<double> &values, bool transpose = false);

    /**
     * Constructs the Vector from the contents of the initializer list @a list.
     * @param list The initializer list.
     */
    Vector(const std::initializer_list<double> &list);

    /** Default copy constructor */
    Vector(const Vector &other) = default;

    /** Default move constructor */
    Vector(Vector &&other) = default;

    /** Default destructor */
    ~Vector() = default;

    /** Default copy assignment operator */
    Vector &operator=(const Vector &other) = default;

    /** Default move assignment operator */
    Vector &operator=(Vector &&other) = default;

    /**
     * @return dimension of vector
     */
    inline count getDimension() const { return values.size(); }

    /**
     * A transposed vector is a row vector.
     * @return True, if this vector is transposed, otherwise false.
     */
    bool isTransposed() const;

    /**
     * @return Transposed copy of this vector.
     */
    Vector transpose() const;

    /**
     * Calculates and returns the Euclidean length of this vector
     * @return The Euclidean length of this vector.
     */
    double length() const;

    /**
     * Calculates and returns the arithmetic mean of this vector
     * @return The arithmetic mean of this vector.
     */
    double mean() const;

    /**
     * Returns a reference to the element at index @a idx without checking the range of this vector.
     * @param idx The index of the element.
     * @return Reference to the element at index @a idx.
     */
    inline double &operator[](index idx) {
        assert(idx < values.size());
        return values[idx];
    }

    /**
     * Assigns the given value to all the elements in the vector.
     *
     * @param val The value to be assigned.
     */
    void fill(double val) { std::fill(values.begin(), values.end(), val); }

    /**
     * Returns the element at index @a idx without checking the range of this vector.
     * @a idx The index of the element.
     * @return Constant reference to the element at index @a idx.
     */
    inline double operator[](index idx) const {
        assert(idx < values.size());
        return values[idx];
    }

    /**
     * Returns a reference to the element at index @a idx. If @a idx is not a valid index an
     * exception is thrown.
     * @param idx The index of the element.
     * @return Reference to the element at index @a idx.
     */
    double &at(index idx) {
        if (idx >= values.size()) {
            throw std::runtime_error("index out of range");
        } else {
            return values[idx];
        }
    }

    /**
     * Compares this vector and @a other element-wise.
     * @return True, if this vector is element-wise equal to @a other, otherwise false.
     */
    bool operator==(const Vector &other) const;

    /**
     * Compares this vector and @a other element-wise.
     * @return True, if this vector is element-wise unequal to @a other, otherwise false.
     */
    bool operator!=(const Vector &other) const;

    /**
     * Computes the outer product of @a v1 and @a v2.
     * @param v1 First Vector.
     * @param v2 Second Vector.
     * @return The resulting matrix from the outer product.
     */
    template <class Matrix = DynamicMatrix>
    static Matrix outerProduct(const Vector &v1, const Vector &v2);

    /**
     * Computes the inner product (dot product) of the vectors @a v1 and @a v2.
     * @return The result of the inner product.
     */
    static double innerProduct(const Vector &v1, const Vector &v2);

    /**
     * Computes the inner product (dot product) of this vector and @a other.
     * @return The result of the inner product.
     */
    double operator*(const Vector &other) const;

    /**
     * Multiplies this vector with @a matrix and returns the result.
     * @return The result of multiplying this vector with @a matrix.
     */
    template <typename Matrix = DynamicMatrix>
    Vector operator*(const Matrix &matrix) const;

    /**
     * Multiplies this vector with a scalar specified in @a scalar and returns the result in a new
     * vector.
     * @return The result of multiplying this vector with @a scalar.
     */
    Vector operator*(double scalar) const;

    /*
     * Multiplies this vector with a scalar specified in @a scalar.
     * @return Reference to this vector.
     */
    Vector &operator*=(double scalar);

    /**
     * Divides this vector by a divisor specified in @a divisor and returns the result in a new
     * vector.
     * @return The result of dividing this vector by @a divisor.
     */
    Vector operator/(double divisor) const;

    /**
     * Divides this vector by a divisor specified in @a divisor.
     * @return Reference to this vector.
     */
    Vector &operator/=(double divisor);

    /**
     * Adds this vector to @a other and returns the result.
     * Note that the dimensions of the vectors have to be the same.
     * @return The sum of this vector and @a other.
     */
    Vector operator+(const Vector &other) const;

    /**
     * Adds @a value to each element of this vector and returns the result.
     */
    Vector operator+(double value) const;

    /**
     * Adds @a other to this vector.
     * Note that the dimensions of the vectors have to be the same.
     * @return Reference to this vector.
     */
    Vector &operator+=(const Vector &other);

    /**
     * Adds @a value to each element of this vector.
     */
    Vector &operator+=(double value);

    /**
     * Subtracts @a other from this vector and returns the result.
     * Note that the dimensions of the vectors have to be the same.
     * @return The difference of this vector and @a other.
     *
     */
    Vector operator-(const Vector &other) const;

    /**
     * Subtracts @a value from each element of this vector and returns the result.
     */
    Vector operator-(double value) const;

    /**
     * Subtracts @a other from this vector.
     * Note that the dimensions of the vectors have to be the same.
     * @return Reference to this vector.
     */
    Vector &operator-=(const Vector &other);

    /**
     * Subtracts @a value from each element of this vector.
     */
    Vector &operator-=(double value);

    /**
     * Applies the unary function @a unaryElementFunction to each value in the Vector. Note that it
     * must hold that the function applied to the zero element of this matrix returns the zero
     * element.
     * @param unaryElementFunction
     */
    template <typename F>
    void apply(F unaryElementFunction);

    /**
     * Iterate over all elements of the vector and call handler (lambda closure).
     */
    template <typename L>
    void forElements(L handle);

    /**
     * Iterate over all elements of the vector and call handler (lambda closure).
     */
    template <typename L>
    void forElements(L handle) const;

    /**
     * Iterate in parallel over all elements of the vector and call handler (lambda closure).
     *
     */
    template <typename L>
    void parallelForElements(L handle);

    /**
     * Iterate in parallel over all elements of the vector and call handler (lambda closure).
     */
    template <typename L>
    void parallelForElements(L handle) const;
};

/**
 * Multiplies the vector @a v with a scalar specified in @a scalar and returns the result.
 * @return The result of multiplying this vector with @a scalar.
 */
inline Vector operator*(double scalar, const Vector &v) {
    return v.operator*(scalar);
}

template <class Matrix>
Matrix Vector::outerProduct(const Vector &v1, const Vector &v2) {
    std::vector<Triplet> triplets;

    for (index i = 0; i < v1.getDimension(); ++i) {
        for (index j = 0; j < v2.getDimension(); ++j) {
            double result = v1[i] * v2[j];
            if (fabs(result) >= FLOAT_EPSILON) {
                triplets.push_back({i, j, result});
            }
        }
    }

    return Matrix(v1.getDimension(), v2.getDimension(), triplets);
}

template <class Matrix>
Vector Vector::operator*(const Matrix &matrix) const {
    assert(isTransposed());                          // vector must be of the form 1xn
    assert(getDimension() == matrix.numberOfRows()); // dimensions of vector and matrix must match

    Vector result(matrix.numberOfColumns(), 0.0, true);
#pragma omp parallel for
    for (omp_index k = 0; k < static_cast<omp_index>(matrix.numberOfColumns()); ++k) {
        Vector column = matrix.column(k);
        result[k] = (*this) * column;
    }

    return result;
}

template <typename F>
void Vector::apply(F unaryElementFunction) {
#pragma omp parallel for
    for (omp_index i = 0; i < static_cast<omp_index>(getDimension()); ++i) {
        values[i] = unaryElementFunction(values[i]);
    }
}

template <typename L>
inline void Vector::forElements(L handle) {
    for (uint64_t i = 0; i < getDimension(); ++i) {
        handle(values[i]);
    }
}

template <typename L>
inline void Vector::forElements(L handle) const {
    for (uint64_t i = 0; i < getDimension(); ++i) {
        handle(values[i]);
    }
}

template <typename L>
inline void Vector::parallelForElements(L handle) {
#pragma omp parallel for
    for (omp_index i = 0; i < static_cast<omp_index>(getDimension()); ++i) {
        handle(i, values[i]);
    }
}

template <typename L>
inline void Vector::parallelForElements(L handle) const {
#pragma omp parallel for
    for (omp_index i = 0; i < static_cast<omp_index>(getDimension()); ++i) {
        handle(i, values[i]);
    }
}

} /* namespace NetworKit */
#endif // NETWORKIT_ALGEBRAIC_VECTOR_HPP_
