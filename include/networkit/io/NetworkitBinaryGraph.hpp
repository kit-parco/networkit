/*
 * NetworkitBinaryGraph.hpp
 *
 *      Author: Charmaine Ndolo <charmaine.ndolo@hu-berlin.de>
 */

#ifndef NETWORKIT_IO_NETWORKIT_BINARY_GRAPH_HPP_
#define NETWORKIT_IO_NETWORKIT_BINARY_GRAPH_HPP_

#include <cstdint>
#include <cstddef>

#include <tlx/math/clz.hpp>
#include <tlx/math/ffs.hpp>
#include <tlx/define/likely.hpp>

namespace NetworKit {
namespace nkbg {

struct Header {
    char magic[8];
    uint64_t checksum;
    uint64_t features;
    uint64_t nodes;
    uint64_t chunks;
    uint64_t offsetBaseData;
    uint64_t offsetAdjLists;
    uint64_t offsetAdjTranspose;
    uint64_t offsetWeightLists;
    uint64_t offsetWeightTranspose;	
    uint64_t offsetAdjIdLists;
    uint64_t offsetAdjIdTranspose;
};

enum WEIGHT_FORMAT {
    NONE = 0,
    VARINT = 1,
    SIGNED_VARINT = 2,
    DOUBLE = 3,
    FLOAT = 4
};

static constexpr uint8_t DELETED_BIT = 0x1; // bit 0
static constexpr uint64_t DIR_MASK = 0x1; // bit 0
static constexpr uint64_t WGHT_MASK = 0xE; //bit 1-3
static constexpr uint64_t WGHT_SHIFT = 0x1;
static constexpr uint64_t INDEX_MASK = 0x10; //bit 4
static constexpr uint64_t INDEX_SHIFT = 0x4;

/**
 * Serializes value into a buffer and returns the number of bytes written.
 * The space required to store an integer x ranges between 1 and 9 bytes.
 *
 * The first byte (buffer[0]) encodes the number of following data bytes
 * in the position y in [0, 7] of the least significant bit that is asserted.
 * If the header byte is zero (i.e., it contains no asserted bits), exactly
 * y=8 data bytes follow.
 * Possibly remaining bits in the header are used to store the least
 * significant bits of x. The subsequent data bytes store the remainder of x
 * in little endian, i.e. * the eight least significant bits (not stored in
 * the header) are stored in buffer[1], the next higher eight bits in
 * buffer[2] and so on.
 *
 * Example:
 * x = 0b0GFE'DCBA with G=1 requires 7 bits of storage. Hence it fits
 * completely into the header which is
 * buffer[0] = 0bGFED'CBA|1 where bit 0 encodes that no
 * data bytes follow while bits 1 to 7 store x.
 *
 * x = 0bHGFE'DCBA with H=1 requires 8 bits of storage. Hence an
 * additional data byte is required:
 *
 * buffer[0] = 0bFEDC'BA|10 where bits 0 and 1 encode that one data byte will follow
 * buffer[1] = 0b0000'00|HG where bits 0 and 1 encode the remaining two bits
 * of x, and bits 2 to 7 are filled as zeros
 */

inline size_t varIntEncode(uint64_t value, uint8_t *buffer) noexcept {
    if (!value) {
        buffer[0] = 1;
        return 1;
    }

    int dataBytes;

    if (TLX_UNLIKELY(value >= (1llu << 56))) {
        // We have more than 56 bits of data, i.e. we need 8 data bytes.
        buffer[0] = 0;
        dataBytes = 8;

    } else {
        const auto bits = 64 - tlx::clz(value);
        dataBytes = static_cast<int>((bits - 1) / 7);

        // number of bytes is encode by the number of trailing zeros
        buffer[0] = static_cast<uint8_t>(1 << dataBytes);

        // put the 8 - (dataBytes - 1) least significant bits into the remainder of the buffer byte
        buffer[0] |= static_cast<uint8_t>(value << (dataBytes + 1));
        value >>= 7 - dataBytes;
    }

    for (int i = 0; i < dataBytes; i++) {
        buffer[i + 1] = static_cast<uint8_t>(value & 0xFF);
        value >>= 8;
    }

    return dataBytes + 1;
}

/// Converts a serialized varint into an uint64_t value and returns the number of bytes consumed.
inline size_t varIntDecode(const uint8_t *data, uint64_t &value) noexcept {
    int n = 8;
    int bitsRecovered = 0;
    uint64_t decoded = 0;

    if (data[0]) {
        n = static_cast<int>(tlx::ffs(data[0]) - 1);
        decoded = data[0] >> (n + 1);
        bitsRecovered = 7 - n;
    }

    for (int i = 0; i < n; i++) {
        decoded |= static_cast<uint64_t>(data[i + 1]) << bitsRecovered;
        bitsRecovered += 8;
    }

    value = decoded;
    return n + 1;
}

/**
 * Encodes a singed value, s.t. the sign bit is in the LSB bit, while the absolute value is
 * stored in the remaining upper 63 bits. This allows more effective compression using varint.
 */
inline uint64_t zigzagEncode(int64_t value) noexcept {
    return (static_cast<uint64_t>(value) << 1) ^ (int64_t(-1) * (value < 0));
}

//! Reverses zigzagEncode.
inline int64_t zigzagDecode(uint64_t value) noexcept {
    return (value >> 1) ^ (-(value & 1));
}

} // namespace nkbg
} // namespace NetworKit

#endif // NETWORKIT_IO_NETWORKIT_BINARY_GRAPH_HPP_
