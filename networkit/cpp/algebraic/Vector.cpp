// no-networkit-format
/*
 * Vector.cpp
 *
 *  Created on: 12.03.2014
 *      Author: Michael Wegner (michael.wegner@student.kit.edu)
 */

#include <networkit/algebraic/Vector.hpp>

#include <networkit/algebraic/DynamicMatrix.hpp>

namespace NetworKit {

Vector::Vector() : values(0), transposed(false) {}

Vector::Vector(const count dimension, const double initialValue, const bool transpose) : values(dimension, initialValue), transposed(transpose) {}

Vector::Vector(const std::vector<double> &values, const bool transpose) : values(values), transposed(transpose) {
}

Vector::Vector(const std::initializer_list<double> &list) : values(list), transposed(false) {
}

bool Vector::isTransposed() const {
    return transposed;
}

Vector Vector::transpose() const {
    Vector v(*this);
    v.transposed = !this->transposed;
    return v;
}

double Vector::length() const {
    return std::sqrt(innerProduct(*this, *this));
}

double Vector::mean() const {
    double sum = 0.0;
    this->forElements([&](double value){
        sum += value;
    });

    return sum / (double) this->getDimension();
}

bool Vector::operator==(const Vector &other) const {
    if (getDimension() != other.getDimension() || isTransposed() != other.isTransposed()) return false;

    for (index i = 0; i < getDimension(); i++) {
        if (values[i] != other[i]) return false;
    }

    return true;
}

bool Vector::operator!=(const Vector &other) const {
    return !(*this == other);
}

double Vector::innerProduct(const Vector &v1, const Vector &v2) {
    assert(v1.getDimension() == v2.getDimension());
    double scalar = 0.0;
    for (index i = 0; i < v1.getDimension(); ++i) {
        scalar += v1[i] * v2[i];
    }

    return scalar;
}

double Vector::operator*(const Vector &other) const {
    assert(isTransposed() && !other.isTransposed()); // vectors must be transposed correctly for inner product
    assert(getDimension() == other.getDimension()); // dimensions of vectors must match

    return innerProduct(*this, other);
}

Vector Vector::operator*(double scalar) const {
    return Vector(*this) *= scalar;
}

Vector& Vector::operator*=(double scalar) {
#pragma omp parallel for
    for (omp_index i = 0; i < static_cast<omp_index>(getDimension()); i++) {
        values[i] *= scalar;
    }

    return *this;
}

Vector Vector::operator/(double divisor) const {
    return Vector(*this) /= divisor;
}

Vector& Vector::operator/=(double divisor) {
    return *this *= 1 / divisor;
}

Vector Vector::operator+(const Vector &other) const {
    return Vector(*this) += other;
}

Vector Vector::operator+(double value) const {
    return Vector(*this) += value;
}

Vector& Vector::operator+=(const Vector &other) {
    assert(isTransposed() == other.isTransposed()); // vectors must be transposed correctly
    assert(getDimension() == other.getDimension()); // dimensions of vectors must match

#pragma omp parallel for
    for (omp_index i = 0; i < static_cast<omp_index>(getDimension()); i++) {
        values[i] += other[i];
    }

    return *this;
}

Vector& Vector::operator+=(const double value) {
#pragma omp parallel for
    for (omp_index i = 0; i < static_cast<omp_index>(getDimension()); ++i) {
        values[i] += value;
    }

    return *this;
}

Vector Vector::operator-(const Vector &other) const {
    return Vector(*this) -= other;
}

Vector Vector::operator-(const double value) const {
    return Vector(*this) -= value;
}

Vector& Vector::operator-=(const Vector &other) {
    assert(isTransposed() == other.isTransposed()); // vectors must be transposed correctly
    assert(getDimension() == other.getDimension()); // dimensions of vectors must match

#pragma omp parallel for
    for (omp_index i = 0; i < static_cast<omp_index>(getDimension()); i++) {
        values[i] -= other[i];
    }

    return *this;
}

Vector& Vector::operator-=(const double value) {
#pragma omp parallel for
    for (omp_index i = 0; i < static_cast<omp_index>(getDimension()); ++i) {
        values[i] -= value;
    }

    return *this;
}


} /* namespace NetworKit */

