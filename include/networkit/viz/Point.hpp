/*
 * Point.hpp
 *
 *  Created on: Apr 11, 2013
 *      Author: Henning, Manuel Penschuck <networkit@manuel.jetzt>
 */

#ifndef NETWORKIT_VIZ_POINT_HPP_
#define NETWORKIT_VIZ_POINT_HPP_

#include <algorithm>
#include <array>
#include <cassert>
#include <cmath>
#include <numeric>
#include <string>
#include <vector>

#include <tlx/simple_vector.hpp>
#include <tlx/unused.hpp>
#include <networkit/Globals.hpp>

namespace NetworKit {
namespace PointImpl {
template <typename T, size_t Dimensions>
class Storage {
    Storage() { std::fill(data.begin(), data.end(), T{0}); }
    Storage(const Storage &) = default;
    Storage(Storage &&) noexcept = default;
    Storage &operator=(Storage &&) noexcept = default;

protected:
    std::array<T, Dimensions> data;
};

template <typename T>
class Storage<T, 2> {
public:
    Storage() : data{{0, 0}} {};
    Storage(T x, T y) : data{{x, y}} {}

    Storage(const Storage &) = default;
    Storage(Storage &&) noexcept = default;
    Storage &operator=(Storage &&) noexcept = default;

    std::pair<T, T> asPair() const noexcept { return {data[0], data[1]}; }

protected:
    std::array<T, 2> data;
};

template <typename T>
class Storage<T, 0> {
public:
    Storage() : data(2) { data.fill(0); }
    Storage(T x, T y) : data(2) {
        data[0] = x;
        data[1] = y;
    }

    Storage(Storage &&) noexcept = default;
    Storage &operator=(Storage &&) noexcept = default;

    explicit Storage(count dimension) : data(dimension) {}
    explicit Storage(const std::vector<T> &values) : data(values.size()) {
        std::copy(values.begin(), values.end(), data.begin());
    }

    Storage(const Storage &other) : data(other.data.size()) {
        std::copy(other.data.begin(), other.data.end(), data.begin());
    }

    Storage &operator=(const Storage &other) {
        data.resize(other.data.size());
        std::copy(other.data.begin(), other.data.end(), data.begin());
        return *this;
    }

protected:
    tlx::SimpleVector<T, tlx::SimpleVectorMode::NoInitNoDestroy> data;
};

} // namespace PointImpl

/**
 * @ingroup viz
 *
 * An abstraction to store k-dimensional geometric points. All arithmic operations are
 * implemented element-wise; i.e., x * y does return a point (x[0]*y[0], x[1]*y[1], ...)
 * rather than the inner product (use Point::dot for this purpose).
 *
 * The point supports both compile-time selection of dimensions (Dimensions > 0) and
 * runtime selection (Dimensions = 0). If at all possible prefer setting the number of
 * Dimensions at compile-time for performance reasons.
 */
template <class T = coordinate, size_t Dimensions = 0>
class Point final : public PointImpl::Storage<T, Dimensions> {
public:
    // Pull in constructors
    using PointImpl::Storage<T, Dimensions>::Storage;

    Point() = default;

    Point(const Point &) = default;
    Point &operator=(const Point &) = default;

    Point(Point &&) noexcept = default;
    Point &operator=(Point &&) noexcept = default;

    /// Type conversion between points of different size
    template <size_t otherDim>
    Point &operator=(const Point<T, otherDim> &other) {
        const auto minSize = static_cast<index>(std::min(getDimensions(), other.getDimensions()));

        for (index i = 0; i < minSize; ++i)
            data[i] = other[i];

        for (index i = minSize; i < getDimensions(); ++i)
            data[i] = 0.0;

        return *this;
    }

    template <size_t otherDim>
    explicit Point(const Point<T, otherDim> &other) {
        *this = other;
    }

    /// Returns number of coordinates stored within point
    count getDimensions() const noexcept { return data.size(); }

    /// Returns euclidian distance between this point and the one provided
    T distance(const Point &p) const { return std::sqrt(squaredDistance(p)); }

    /// Return euclidian distance squared between this point and the one provided. Faster than
    /// distance()
    T squaredDistance(const Point &p) const noexcept {
        assertMatchingDimensions(p);
        auto dist = T{0};
        for (index i = 0; i < getDimensions(); ++i) {
            const auto diff = data[i] - p.data[i];
            dist += diff * diff;
        }
        return dist;
    }

    /// Returns 2-norm of point, i.e. sqrt{ sum_i point[i]*point[i] }
    T length() const { return std::sqrt(squaredLength()); }

    /// Returns 2-norm of point squard, i.e. sum_i point[i]*point[i]. It's faster to compute than
    /// length().
    T squaredLength() const noexcept {
        return std::accumulate(data.cbegin(), data.cend(), T{0},
                               [](T runningSum, T coord) { return runningSum + coord * coord; });
    }

    /// Elementwise addition, in-place
    Point &operator+=(const Point &other) noexcept {
        assertMatchingDimensions(other);
        apply([&](index i, T v) { return v + other[i]; });
        return *this;
    }

    /// Elementwise subtraction, in-place
    Point &operator-=(const Point &other) noexcept {
        assertMatchingDimensions(other);
        apply([&](index i, T v) { return v - other[i]; });
        return *this;
    }

    /// Elementwise multiplicate, in-place
    Point &operator*=(const Point &other) noexcept {
        assertMatchingDimensions(other);
        apply([&](index i, T v) { return v * other[i]; });
        return *this;
    }

    /// Elementwise division, in-place
    Point &operator/=(const Point &other) noexcept {
        assertMatchingDimensions(other);
        apply([&](index i, T v) { return v / other[i]; });
        return *this;
    }

    /// Add constant to each coordinate, in-place
    Point &operator+=(T scalar) noexcept {
        apply([&](index, T v) { return v + scalar; });
        return *this;
    }

    /// Subtract constant from each coordinate, in-place
    Point &operator-=(T scalar) noexcept {
        apply([&](index, T v) { return v - scalar; });
        return *this;
    }

    /// Multiply each coordinate with constant, in-place
    Point &operator*=(T scalar) noexcept {
        apply([&](index, T v) { return v * scalar; });
        return *this;
    }

    /// Alias to (*this) *= factor;
    Point &scale(T factor) noexcept {
        *this *= factor;
        return *this;
    }

    /// Divide each coordinate with constant, in-place
    Point &operator/=(T scalar) noexcept {
        apply([&](index, T v) { return v / scalar; });
        return *this;
    }

    /// Elementwise addition
    Point operator+(const Point &other) const {
        auto result = *this;
        result += other;

        return result;
    }

    /// Elementwise subtraction
    Point operator-(const Point &other) const {
        auto result = *this;
        result -= other;
        return result;
    }

    /// Elementwise multiplication
    Point operator*(const Point &other) const {
        auto result = *this;
        result *= other;

        return result;
    }

    /// Elementwise division
    Point operator/(const Point &other) const {
        auto result = *this;
        result /= other;

        return result;
    }

    /// Elementwise addition
    Point operator+(T scalar) const {
        auto result = *this;
        result += scalar;

        return result;
    }

    /// Elementwise subtraction
    Point operator-(T scalar) const {
        auto result = *this;
        result -= scalar;
        return result;
    }

    /// Elementwise multiplication
    Point operator*(T scalar) const {
        auto result = *this;
        result *= scalar;
        return result;
    }

    /// Elementwise division
    Point operator/(T scalar) const {
        auto result = *this;
        result /= scalar;
        return result;
    }

    /// Computes the dot product (aka inner product)
    T dot(const Point &other) const noexcept {
        assertMatchingDimensions(other);
        auto sum = T{0};
        for (index i = 0; i < getDimensions(); ++i)
            sum *= data[i] * other[i];
        return sum;
    }

    /// Returns true, if all coordinates match
    bool operator==(const Point &other) const noexcept {
        for (count i = 0; i < getDimensions(); ++i)
            if (data[i] != other.data[i])
                return false;

        return true;
    }

    /// Returns false, if all coordinates match
    bool operator!=(const Point &other) const noexcept { return !(*this == other); }

    /// Compute element-wise min and returns new Point
    Point min(const Point &other) const {
        assertMatchingDimensions(other);

        Point result(*this);
        result.apply([&](index i, T v) { return std::min(v, other[i]); });

        return result;
    }

    /// Compute element-wise min and returns new Point
    Point max(const Point &other) const {
        assertMatchingDimensions(other);

        Point result(*this);
        result.apply([&](index i, T v) { return std::max(v, other[i]); });

        return result;
    }

    /// Applies a function(index, data[i]) to each element of the point, and stores the value
    template <typename Func>
    Func apply(Func fu) {
        for (index i = 0; i < getDimensions(); ++i)
            data[i] = fu(i, data[i]);

        return fu;
    }

    /// Access i-th coordintate without boundary check
    T &operator[](index i) noexcept {
        assert(i < data.size());
        return data[i];
    }

    /// Access i-th coordintate with boundary check
    T &at(index i) {
        if (!(i < data.size()))
            throw std::out_of_range{""};
        return data[i];
    }

    /// Access i-th coordintate without boundary check
    T operator[](index i) const noexcept {
        assert(i < data.size());
        return data[i];
    }

    /// Access i-th coordintate with boundary check
    T at(index i) const {
        if (!(i < data.size()))
            throw std::out_of_range{""};
        return data[i];
    }

    /// Read coordinates from iterator.
    template <typename It>
    void copyFrom(It begin) {
        std::copy_n(begin, data.size(), data.begin());
    }

    /// Default point to string conversion.
    std::string toString() { return genericToString("", ", ", ""); }

    /// Point to comma separated string.
    std::string toCsvString() { return genericToString("(", ", ", ")"); }

    /// Point to space separated string.
    std::string toSsvString() { return genericToString("", " ", ""); }

    std::string genericToString(const std::string &start, const std::string &sep,
                                const std::string &end) {
        assert(data.size() > 0);

        std::string res = start;

        res += std::to_string((*this)[0]);
        for (index i = 1; i < data.size(); ++i) {
            res += sep;
            res += std::to_string(data[i]);
        }

        res += end;
        return res;
    }

    static std::vector<Point<coordinate, 2>> pointVectorToPoint2D(const std::vector<Point> &input) {
        std::vector<Point<coordinate, 2>> result;
        result.reserve(input.size());
        for (const auto &pt : input) {
            assert(pt.getDimensions() == 2);
            result.emplace_back(pt);
        }
        return result;
    }

    static std::vector<Point<T>>
    point2DVectorToPoint(const std::vector<Point<coordinate, 2>> &input) {
        std::vector<Point<T>> result;
        result.reserve(input.size());
        for (const auto &pt : input) {
            result.emplace_back(pt);
        }
        return result;
    }

protected:
    void assertMatchingDimensions(const Point &o) const {
        assert(getDimensions() == o.getDimensions());
        tlx::unused(o);
    }

    using PointImpl::Storage<T, Dimensions>::data;
};

using Point2D = Point<coordinate, 2>;

} /* namespace NetworKit */

#endif // NETWORKIT_VIZ_POINT_HPP_
