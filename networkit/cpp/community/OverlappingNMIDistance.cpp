
#include <cmath>
#include <iomanip>
#include <numeric>

#include <networkit/auxiliary/NumericTools.hpp>
#include <networkit/auxiliary/SignalHandling.hpp>
#include <networkit/community/OverlappingNMIDistance.hpp>

namespace NetworKit {

OverlappingNMIDistance::SizesAndIntersections
OverlappingNMIDistance::calculateClusterAndIntersectionSizes(const Graph &graph, const Cover &X,
                                                             const Cover &Y) {
    SizesAndIntersections result;
    result.sizesX.resize(X.upperBound());
    result.sizesY.resize(Y.upperBound());

    // nodeRange() ignores deleted nodes.
    for (node u : graph.nodeRange()) {
        for (auto i : X[u]) {
            for (auto j : Y[u]) {
                ++result.intersectionSizes[{i, j}];
            }

            ++result.sizesX[i];
        }
        for (auto j : Y[u]) {
            ++result.sizesY[j];
        }
    }

    return result;
}

double OverlappingNMIDistance::h(count w, count n) {
    auto v = static_cast<double>(w);
    auto value = v > 0.0 ? -v * std::log2(v / n) : 0.0;
    assert(!std::isnan(value));
    return value;
}

double OverlappingNMIDistance::entropy(count size, count n) {
    assert(size <= n);
    auto value = h(size, n) + h(n - size, n);
    assert(value >= 0.0);
    return value;
}

double OverlappingNMIDistance::entropy(const std::vector<count> &sizesX, count n) {
    double entropyX = 0.0;
    for (auto sizeXi : sizesX) {
        if (sizeXi == 0)
            continue;
        entropyX += entropy(sizeXi, n);
    }
    assert(!std::isnan(entropyX));
    assert(entropyX >= 0.0);
    return entropyX;
}

double OverlappingNMIDistance::adjustedConditionalEntropy(count sizeXi, count sizeYj,
                                                          count intersectionSize, count n) {
    count a = n + intersectionSize - sizeXi - sizeYj; // X_{i,u} = 0 and Y_{j,u} = 0
    count b = sizeYj - intersectionSize;              // X_{i,u} = 0 and Y_{j,u} = 1
    count c = sizeXi - intersectionSize;              // X_{i,u} = 1 and Y_{j,u} = 0
    count d = intersectionSize;                       // X_{i,u} = 1 and Y_{j,u} = 1

    auto entropyXiYj = h(a, n) + h(b, n) + h(c, n) + h(d, n);
    auto entropyYj = entropy(sizeYj, n);
    auto entropyXi = entropy(sizeXi, n);

    if (h(a, n) + h(d, n) >= h(b, n) + h(c, n)) {
        return entropyXiYj - entropyYj;
    } else {
        return entropyXi;
    }
}

double OverlappingNMIDistance::conditionalEntropy(
    const std::vector<count> &sizesX, const std::vector<count> &sizesY,
    const std::unordered_map<std::pair<index, index>, count, Aux::PairHash> &intersectionSizes,
    bool invertPairIndices, count n) {

    // Choice to initialize H(X_i|Y) with H(X_i), same as reference implementation
    // https://github.com/aaronmcdaid/Overlapping-NMI/blob/master/onmi.cpp#L232
    // This is different to the mathematical definition.
    // For example: X = {{1..80}}, Y = {{81..83}}, n = 100
    // H^*(X_1|Y_1) ~= 64.95 is lower than H(X_i) ~= 72.19, but H^*(X_1|Y_1) is not calculated,
    // because X_1 and Y_1 do not intersect.
    // NOTE: Handle clusters which do not intersect any cluster of the other cover
    //       Check if a cluster X_i does not intersect with any cluster Y_j from the other cover.
    //       For Y_j in increasing size calculate H^*(X_i|Y_j) to find the minimum.
    //       Stop as soon as h(n - |X_i| - |Y_j|, n) >= h(|Y_j|, n) + h(|X_i|, n) no longer holds.

    // Stores H(X_i|Y) = \min_{Y_j \in Y} H^*(X_i|Y_j) for each i
    std::vector<double> conditionalEntropiesXiGivenY(sizesX.size(), 0.0);
    for (index i = 0; i < sizesX.size(); ++i) {
        if (sizesX[i] == 0)
            continue;
        conditionalEntropiesXiGivenY[i] = entropy(sizesX[i], n);
    }

    for (const auto &entry : intersectionSizes) {
        auto i = entry.first.first;
        auto j = entry.first.second;
        if (invertPairIndices)
            std::swap(i, j);
        auto intersectionSize = entry.second;

        conditionalEntropiesXiGivenY[i] =
            std::min(conditionalEntropiesXiGivenY[i],
                     adjustedConditionalEntropy(sizesX[i], sizesY[j], intersectionSize, n));
    }

    auto conditionalEntropyXGivenY = std::accumulate(conditionalEntropiesXiGivenY.begin(),
                                                     conditionalEntropiesXiGivenY.end(), 0.0);

    assert(!std::isnan(conditionalEntropyXGivenY));
    return conditionalEntropyXGivenY;
}

void OverlappingNMIDistance::clampBelow(double &value, double lowerBound, const char *const format,
                                        int printPrecision) {
    if (value < lowerBound) {
        if (!Aux::NumericTools::ge(value, lowerBound, Aux::NumericTools::acceptableError)) {
            std::stringstream ss;
            ss << std::fixed << std::setprecision(printPrecision) << value;
            const auto s = ss.str();
            ERRORF(format, s.c_str());
        }
        value = lowerBound;
    }
}

void OverlappingNMIDistance::clampAbove(double &value, double upperBound, const char *const format,
                                        int printPrecision) {
    if (value > upperBound) {
        if (!Aux::NumericTools::le(value, upperBound, Aux::NumericTools::acceptableError)) {
            std::stringstream ss;
            ss << std::fixed << std::setprecision(printPrecision) << value;
            const auto s = ss.str();
            ERRORF(format, s.c_str());
        }
        value = upperBound;
    }
}

double OverlappingNMIDistance::normalize(OverlappingNMIDistance::Normalization normalization,
                                         double mutualInformation, double entropyX,
                                         double entropyY) {
    // Numerical errors less than Aux::NumericTools::acceptableError are corrected silently.
    // Larger errors are also corrected, but they are logged as an ERROR.
    clampBelow(entropyX, 0.0, "Set entropyX lower than 0.0 to 0.0: %s");
    clampBelow(entropyY, 0.0, "Set entropyY lower than 0.0 to 0.0: %s");

    // Assumption: Empty covers were already dealt with.
    // At this point, H(X) = 0 implies that X has only clusters with all nodes.
    if (entropyX == 0.0 && entropyY == 0.0) {
        // Duplicate clusters should not exist, but its not explicitly checked.
        // De-duplicating clusters results in two covers that are equal.
        return 1.0;
    } else if ((entropyX == 0.0 || entropyY == 0.0)
               && (normalization == OverlappingNMIDistance::Normalization::MIN
                   || normalization == OverlappingNMIDistance::Normalization::GEOMETRIC_MEAN)) {
        // Prevent division by zero.
        // One of the covers has only clusters with all nodes. The other cover is different.
        return 0.0;
    }

    double nmi;
    switch (normalization) {
    case OverlappingNMIDistance::Normalization::MIN:
        nmi = mutualInformation / std::min(entropyX, entropyY);
        break;
    case OverlappingNMIDistance::Normalization::GEOMETRIC_MEAN:
        nmi = mutualInformation / std::sqrt(entropyX * entropyY);
        break;
    case OverlappingNMIDistance::Normalization::ARITHMETIC_MEAN:
        nmi = 2 * mutualInformation / (entropyX + entropyY);
        break;
    case OverlappingNMIDistance::Normalization::MAX:
        nmi = mutualInformation / std::max(entropyX, entropyY);
        break;
    case OverlappingNMIDistance::Normalization::JOINT_ENTROPY: {
        auto entropyXY = entropyX + entropyY - mutualInformation;
        nmi = mutualInformation / entropyXY;
        break;
    }
    default:
        throw std::logic_error("normalization method is not covered");
    }

    if (std::isnan(nmi)) {
        ERROR("Set nmi ", nmi, " to 0.0");
        nmi = 0.0;
    }
    clampBelow(nmi, 0.0, "Set nmi lower than 0.0 to 0.0: %s");
    clampAbove(nmi, 1.0, "Set nmi larger than 1.0 to 1.0: %s");

    return nmi;
}

double OverlappingNMIDistance::getDissimilarity(const Graph &G, const Partition &zeta,
                                                const Partition &eta) {
    return getDissimilarity(G, Cover(zeta), Cover(eta));
}

double OverlappingNMIDistance::getDissimilarity(const Graph &G, const Cover &zeta,
                                                const Cover &eta) {

    // Assumptions:
    //   All clusters are non-empty. When encountering one (due to non-consecutive ids), they are
    //   ignored.
    //   The clusters of a cover are unique. The code does not check that explicitly. The assumption
    //   is used when dealing with edge cases.

    // Edge cases:
    //   For H(X) = 0, X is either an empty cover or X has only clusters with all nodes.
    //   This leads to numerical problems with the normalization.

    Aux::SignalHandler handler;

    const auto n = G.numberOfNodes();

    if (zeta.numberOfElements() != G.upperNodeIdBound()
        || eta.numberOfElements() != G.upperNodeIdBound()) {
        // numberOfElements() and upperNodeIdBound() both also count deleted nodes.
        throw std::invalid_argument("Covers must have the same number of nodes as the graph.");
    }

    const auto sizes = calculateClusterAndIntersectionSizes(G, zeta, eta);
    const auto &sizesX = sizes.sizesX;
    const auto &sizesY = sizes.sizesY;
    const auto &intersectionSizes = sizes.intersectionSizes;

    // Edge cases for empty covers.
    const auto isEmpty = [](count size) { return size == 0; };
    const auto XisEmpty = std::all_of(sizesX.begin(), sizesX.end(), isEmpty);
    const auto YisEmpty = std::all_of(sizesY.begin(), sizesY.end(), isEmpty);
    if (XisEmpty != YisEmpty) {
        // The covers are different and one is empty.
        return 1.0;
    } else if (XisEmpty && YisEmpty) {
        // Both covers are empty.
        return 0.0;
    }

    handler.assureRunning();

    auto conditionalEntropyXGivenY =
        conditionalEntropy(sizesX, sizesY, intersectionSizes, false, n);
    auto conditionalEntropyYGivenX = conditionalEntropy(sizesY, sizesX, intersectionSizes, true, n);

    auto entropyX = entropy(sizesX, n);
    auto entropyY = entropy(sizesY, n);

    auto mutualInformation =
        0.5 * (entropyX - conditionalEntropyXGivenY + entropyY - conditionalEntropyYGivenX);

    auto nmi = normalize(mNormalization, mutualInformation, entropyX, entropyY);
    auto nmiDistance = 1.0 - nmi;

    assert(0.0 <= nmiDistance && nmiDistance <= 1.0);

    return nmiDistance;
}

} // namespace NetworKit
