// no-networkit-format
/*
 * Point2DWithIndex.h
 *
 *  Created on: 24.07.2014
 *      Author: moritzl
 */

#ifndef NETWORKIT_GEOMETRIC_POINT2_D_WITH_INDEX_HPP_
#define NETWORKIT_GEOMETRIC_POINT2_D_WITH_INDEX_HPP_

#include <cassert>
#include <cinttypes>
#include <cmath>
#include <cstdint>
#include <vector>

#include <networkit/Globals.hpp>

namespace NetworKit {

/**
 * @ingroup viz
 * Points in any dimension of templated type.
 */
template<class T>
class Point2DWithIndex final {
private:
    T x;
    T y;
    index indice;

public:
    Point2DWithIndex() {
        x = 0;
        y = 0;
        indice = 0;
    }
    Point2DWithIndex(T x, T y) {
        this->x = x;
        this->y = y;
        this->indice = 0;
    }
    Point2DWithIndex(T x, T y, index indice) {
        this->x = x;
        this->y = y;
        this->indice = indice;
    }

    count getDimensions() { return 2; }

    T distance(const Point2DWithIndex<T>& p) const;
    T squaredDistance(const Point2DWithIndex<T>& p) const;

    Point2DWithIndex& operator+=(const Point2DWithIndex<T>& p);
    Point2DWithIndex& operator-=(const Point2DWithIndex<T>& p);
    Point2DWithIndex& scale(T factor);

    Point2DWithIndex operator-(const Point2DWithIndex<T>& other);
    Point2DWithIndex operator+(const Point2DWithIndex<T>& other);

    T length() const;
    T squaredLength() const;

    T& operator[](index i);
    T getX() const;
    T getY() const;
    index getIndex() const;
};

template<class T>
T Point2DWithIndex<T>::length() const {
    return sqrt(x*x+y*y);
}

template<class T>
T Point2DWithIndex<T>::squaredLength() const {
    return x*x + y*y;
}

template<class T>
T Point2DWithIndex<T>::squaredDistance(const Point2DWithIndex<T>& p) const {
    T diffx = p.x - x;
    T diffy = p.y - y;
    return diffx*diffx + diffy*diffy;
}

template<class T>
T Point2DWithIndex<T>::distance(const Point2DWithIndex<T>& p) const {
    return sqrt(squaredDistance(p));
}

template<class T>
Point2DWithIndex<T>& Point2DWithIndex<T>::operator+=(const Point2DWithIndex<T>& p) {
    this->x += p.x;
    this->y += p.y;
    return *this;
}

template<class T>
Point2DWithIndex<T>& Point2DWithIndex<T>::operator-=(const Point2DWithIndex<T>& p) {
    this->x -= p.x;
    this->y -= p.y;
    return *this;
}

template<class T>
Point2DWithIndex<T> Point2DWithIndex<T>::operator-(const Point2DWithIndex<T>& other) {
    return Point2DWithIndex(x - other.x, y - other.y);
}

template<class T>
Point2DWithIndex<T> Point2DWithIndex<T>::operator+(const Point2DWithIndex<T>& other) {
    return Point2DWithIndex(x + other.x, y + other.y);
}


template<class T>
Point2DWithIndex<T>& Point2DWithIndex<T>::scale(T factor) {
    x *= factor;
    y *= factor;
    return *this;
}

template<class T>
inline T& Point2DWithIndex<T>::operator [](index i) {
    assert(i < 2);
    return i > 0 ? y : x;
}

template<class T>
inline T Point2DWithIndex<T>::getX() const {
    return x;
}

template<class T>
inline T Point2DWithIndex<T>::getY() const {
    return y;
}

template<class T>
inline index Point2DWithIndex<T>::getIndex() const {
    return indice;
}

} /* namespace NetworKit */
#endif // NETWORKIT_GEOMETRIC_POINT2_D_WITH_INDEX_HPP_
