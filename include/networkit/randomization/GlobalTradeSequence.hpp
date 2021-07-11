/*
 * GlobalTradeSequence.hpp
 *
 * This header file is deprecated and will eventually be removed.
 * Do not included it directly.
 *
 *  Created on: 23.05.2018
 *      Author: Manuel Penschuck <networkit@manuel.jetzt>
 */

#ifndef NETWORKIT_RANDOMIZATION_GLOBAL_TRADE_SEQUENCE_HPP_
#define NETWORKIT_RANDOMIZATION_GLOBAL_TRADE_SEQUENCE_HPP_

#ifndef NETWORKIT_PRIVATE_RANDOMIZATION_GLOBAL_TRADE_SEQUENCE_HPP_
#if defined(__clang__) || defined(__GNUG__)
#warning "This header file is deprecated. Do not included it directly."
#elif defined(_MSC_VER)
#pragma message("This header file is deprecated. Do not included it directly.")
#endif
#endif // NETWORKIT_PRIVATE_RANDOMIZATION_GLOBAL_TRADE_SEQUENCE_HPP_

#include <cassert>
#include <cmath>
#include <iostream>
#include <random>
#include <tuple>
#include <type_traits>

#include <networkit/Globals.hpp>
#include <networkit/auxiliary/Log.hpp>

namespace NetworKit {
namespace CurveballDetails {

template <typename T>
class FixedLinearCongruentialMap; // forward declaration (implementation below)

/**
 * Computes a bijection f:[p]->[p], x -> (a*x+b) mod p where p is chosen
 * as the smallest p >= n and p prime. a and b are drawn unif at random.
 * @tparam T
 */
template <typename T>
class LinearCongruentialMap {
    static_assert(!std::is_signed<T>::value, "LinearCongruentialMap requires unsigned types");

    using signed_value_type = typename std::make_signed<T>::type;
    using signed_tuple = std::tuple<signed_value_type, signed_value_type, signed_value_type>;

    friend FixedLinearCongruentialMap<T>;

public:
    using value_type = T;

    LinearCongruentialMap() {}

    explicit LinearCongruentialMap(value_type n, unsigned seed = 0) : n(n), p(computeNextPrime(n)) {
        std::mt19937_64 prng(seed * n + n);
        sampleParameters(prng);
    }

    LinearCongruentialMap(value_type n, std::mt19937_64 &prng) : n(n), p(computeNextPrime(n)) {
        sampleParameters(prng);
    }

    LinearCongruentialMap(value_type n, value_type a, value_type b)
        : n(n), p(computeNextPrime(n)), a(a),
          ainv(static_cast<value_type>((std::get<1>(gcdExtended(a, p)) + p) % p)), b(b) {}

    LinearCongruentialMap(const LinearCongruentialMap &) = default;
    LinearCongruentialMap(LinearCongruentialMap &&) noexcept = default;
    LinearCongruentialMap &operator=(const LinearCongruentialMap &) = default;
    LinearCongruentialMap &operator=(LinearCongruentialMap &&) noexcept = default;

    //! Hashes [n] to [p] with (a*x+b) mod p
    value_type hash(value_type n) const { return (a * n + b) % p; }

    //! Alias to hash(n)
    value_type operator()(value_type n) const { return hash(n); }

    //! invert(hash(x)) == x
    value_type invert(value_type y) const { return ainv * (y + p - b) % p; }

    //! Is true if there exists no x in [n], s.t. h(x) == y,
    //! i.e. it can be used to determine if the map maps to y
    bool isGap(value_type y) const { return invert(y) >= n; }

    //! randomly samples parameters a and b
    void sampleParameters(std::mt19937_64 &prng) {
        const value_type max_a = std::numeric_limits<value_type>::max() / n - 1;
        if (max_a < p) {
            std::cerr << "WARNING: Reduce randomness of hash function to avoid integer "
                         "precision issues\n";
        }

        a = std::uniform_int_distribution<value_type>{1, std::min<value_type>(p, max_a) - 1}(prng);
        ainv = static_cast<value_type>((std::get<1>(gcdExtended(a, p)) + p) % p);
        b = std::uniform_int_distribution<value_type>{0, p - 1}(prng);
    }

    //! Sets parameters a = 1 and b = 0
    void setAsIdentity() {
        a = 1;
        b = 0;
        ainv = 1;
    }

    value_type param_a() const { return a; }
    value_type param_ainv() const { return ainv; }
    value_type param_b() const { return b; }
    value_type param_p() const { return p; }

private:
    value_type n; //< number of elements to be mapped
    value_type p; //< size of the universe chosen to be >= n

    value_type a;    //< multiplicative parameter in hash
    value_type ainv; //< multiplicative invert

    value_type b; //< additive parameter in hash

    value_type computeNextPrime(value_type n) const {
        auto is_prime = [](value_type n) {
            if (n <= 3)
                return true;
            if (0 == n % 2 || 0 == n % 3)
                return false;

            const value_type sqrt = static_cast<value_type>(std::sqrt(n) + 2);
            for (value_type i = 5; i < sqrt; i += 6)
                if (0 == n % i)
                    return false;

            return true;
        };

        while (!is_prime(n))
            n++;

        return n;
    }

    // extended Euclidean algorithm with 1 = gcd(a, p) = a*s + t*p mod p = a*s
    // --> s = 1/a
    static signed_tuple gcdExtended(signed_value_type a, signed_value_type b) noexcept {
        if (a == 0)
            return signed_tuple{b, 0, 1};

        const value_type div = b / a;
        const value_type rem = b % a;

        const auto recursion = gcdExtended(rem, a);
        value_type x = std::get<2>(recursion) - div * std::get<1>(recursion);

        return signed_tuple(std::get<0>(recursion), x, std::get<1>(recursion));
    }
};

/**
 * Computes a bijection f:[p]->[p], x -> (a*x+b) mod p where p is chosen
 * as the smallest p >= n and p prime. a and b are drawn unif at random.
 * @tparam T
 */
template <typename T>
class FixedLinearCongruentialMap {
    static_assert(!std::is_signed<T>::value, "LinearCongruentialMap requires unsigned types");

    using signed_value_type = typename std::make_signed<T>::type;
    using signed_tuple = std::tuple<signed_value_type, signed_value_type, signed_value_type>;

public:
    using value_type = T;

    FixedLinearCongruentialMap() {}

    explicit FixedLinearCongruentialMap(value_type n, unsigned seed = 0) : n(n) {
        std::mt19937_64 prng(seed * n + n);
        sampleParameters(prng);
    }

    FixedLinearCongruentialMap(value_type n, std::mt19937_64 &prng) : n(n) {
        sampleParameters(prng);
    }

    FixedLinearCongruentialMap(value_type n, value_type a, value_type b)
        : n(n), a(a), //
          ainv(static_cast<value_type>(
              (std::get<1>(LinearCongruentialMap<T>::gcdExtended(a, p)) + p) % p)),
          b(b) {}

    FixedLinearCongruentialMap(const FixedLinearCongruentialMap &) = default;
    FixedLinearCongruentialMap(FixedLinearCongruentialMap &&) noexcept = default;
    FixedLinearCongruentialMap &operator=(const FixedLinearCongruentialMap &) = default;
    FixedLinearCongruentialMap &operator=(FixedLinearCongruentialMap &&) noexcept = default;

    //! Hashes [n] to [p] with (a*x+b) mod p
    value_type hash(value_type n) const { return (a * n + b) % p; }

    //! Alias to hash(n)
    value_type operator()(value_type n) const { return hash(n); }

    //! invert(hash(x)) == x
    value_type invert(value_type y) const { return ainv * (y + p - b) % p; }

    //! Is true if there exists no x in [n], s.t. h(x) == y,
    //! i.e. it can be used to determine if the map maps to y
    bool isGap(value_type y) const { return invert(y) >= n; }

    //! randomly samples parameters a and b
    void sampleParameters(std::mt19937_64 &prng) {
        if (n >= p)
            throw std::runtime_error("Support only up to 2147483646 nodes");

        a = std::uniform_int_distribution<value_type>{1, p - 1}(prng);
        ainv = static_cast<value_type>(
            (std::get<1>(LinearCongruentialMap<T>::gcdExtended(a, p)) + p) % p);
        b = std::uniform_int_distribution<value_type>{0, p - 1}(prng);
    }

    //! Sets parameters a = 1 and b = 0
    void setAsIdentity() {
        a = 1;
        b = 0;
        ainv = 1;
    }

    value_type param_a() const { return a; }
    value_type param_ainv() const { return ainv; }
    value_type param_b() const { return b; }
    value_type param_p() const { return p; }

private:
    value_type n; //< number of elements to be mapped
    static constexpr value_type p = 2147483647;

    value_type a;    //< multiplicative parameter in hash
    value_type ainv; //< multiplicative invert

    value_type b; //< additive parameter in hash
};

/**
 * Handles a sequence of num_global_trades hash functions each with its own
 * set of parameters. It prominently exposes the "current" hash function
 * and its successors. Switching the current hash function is possible with
 * using the switch_to_round method.
 *
 * hash, invert  access the "current" hash function
 * hash_next     accesses the successor of "current" which is the identity if
 *               there are no further hash functions
 *
 * @tparam Default constructable, copy-able, requires methods
 *         sample_parameters, set_as_identity, hash and (optinally) invert
 */
template <typename Hash>
class GlobalTradeSequence {
public:
    using value_type = typename Hash::value_type;

    GlobalTradeSequence(node num_nodes, size_t num_global_trades, std::mt19937_64 &prng) {
        hashFunctors.reserve(num_global_trades);
        hashFunctors.push_back(Hash{num_nodes, prng});

        while (hashFunctors.size() < num_global_trades) {
            // we copy the hash function, rather than constructing a new one,
            // to avoid repeated computations, such as the next larger prime
            // number
            hashFunctors.push_back(hashFunctors.back());
            hashFunctors.back().sampleParameters(prng);
        }

        switchToRound(0);
    }

    void switchToRound(size_t round) {
        assert(round < hashFunctors.size());
        current = hashFunctors[round];
        if (round + 1 == hashFunctors.size()) {
            next = current;
            next.setAsIdentity();
        } else {
            next = hashFunctors[round + 1];
        }
    }

    value_type hash(node u) const { return current.hash(u); }

    node invert(value_type u) const { return current.invert(u); }

    value_type hashNext(node u) const { return next.hash(u); }

    size_t numberOfRounds() const { return hashFunctors.size(); }

private:
    std::vector<Hash> hashFunctors;
    Hash current;
    Hash next;
};

} // namespace CurveballDetails
} // namespace NetworKit

#endif // NETWORKIT_RANDOMIZATION_GLOBAL_TRADE_SEQUENCE_HPP_
