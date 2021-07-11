// no-networkit-format
#include <algorithm>
#include <cmath>

#include <networkit/auxiliary/Log.hpp>
#include <networkit/auxiliary/NumericTools.hpp>
#include <networkit/auxiliary/Random.hpp>
#include <networkit/generators/PowerlawDegreeSequence.hpp>

NetworKit::PowerlawDegreeSequence::PowerlawDegreeSequence(NetworKit::count minDeg, NetworKit::count maxDeg, double gamma) :
    minDeg(minDeg), maxDeg(maxDeg), gamma(gamma) {
    if (minDeg > maxDeg) throw std::runtime_error("Error: minDeg must not be larger than maxDeg");
    if (gamma > -1) throw std::runtime_error("Error: gamma must be lower than -1");
}

NetworKit::PowerlawDegreeSequence::PowerlawDegreeSequence(const std::vector< double > &degreeSequence) : minDeg(std::numeric_limits<count>::max()), maxDeg(std::numeric_limits<count>::min()) {
    count sum = 0;
    for (auto &d : degreeSequence) {
        if (d < minDeg) minDeg = d;
        if (d > maxDeg) maxDeg = d;
        sum += d;
    }

    double avg = sum * 1.0 / degreeSequence.size();

    setGammaFromAverageDegree(avg);
}

NetworKit::PowerlawDegreeSequence::PowerlawDegreeSequence(const NetworKit::Graph &g) : minDeg(std::numeric_limits<count>::max()), maxDeg(std::numeric_limits<count>::min()) {
    count sum = 0;
    g.forNodes([&](node u) {
        count d = g.degree(u);
        if (d < minDeg) minDeg = d;
        if (d > maxDeg) maxDeg = d;
        sum += d;
    });

    double avg = sum * 1.0 / g.numberOfNodes();

    setGammaFromAverageDegree(avg);
}



void NetworKit::PowerlawDegreeSequence::setMinimumDegree(NetworKit::count minDeg) {
    this->minDeg = minDeg;
    hasRun = false;
}

void NetworKit::PowerlawDegreeSequence::setGamma(double gamma) {
    this->gamma = gamma;
    hasRun = false;
}

void NetworKit::PowerlawDegreeSequence::setMinimumFromAverageDegree(double avgDeg) {
    double dmin_l = 1.0;
    double dmin_r = maxDeg;
    setMinimumDegree(dmin_l); run();
    double average_l = getExpectedAverageDegree();
    double average_r = maxDeg;

    if (average_l > avgDeg) {
        throw std::runtime_error("The average degree is too low");
    }

    if (average_r < avgDeg) {
        throw std::runtime_error("The average degree must not be higher than the maximum degree");
    }

    while (dmin_l + 1 < dmin_r) {
        setMinimumDegree(static_cast<double>(dmin_r + dmin_l) * 0.5); run();
        double avg = getExpectedAverageDegree();

        TRACE("Trying minDeg ", minDeg, ", this gives average ", avg, ", which should be between ", average_r, " and ", average_l);

        if (avg > avgDeg) {
            average_r = avg;
            dmin_r = minDeg;
        } else {
            average_l = avg;
            dmin_l = minDeg;
        }
    }

    if (avgDeg - average_l < average_r - avgDeg) {
        minDeg = dmin_l;
    } else {
        minDeg = dmin_r;
    }

    hasRun = false;
}

void NetworKit::PowerlawDegreeSequence::setGammaFromAverageDegree(double avgDeg, double minGamma, double maxGamma) {
    double gamma_l = maxGamma;
    double gamma_r = minGamma;
    setGamma(gamma_l); run();
    double average_l = getExpectedAverageDegree();
    setGamma(gamma_r); run();
    double average_r = getExpectedAverageDegree();

    // Note: r is the larger expected average degree!
    if (avgDeg > average_r) {
        setGamma(gamma_r);
        return;
    }

    if (avgDeg < average_l) {
        setGamma(gamma_l);
        return;
    }

    while (gamma_l + 0.001 < gamma_r) {
        setGamma((gamma_r + gamma_l) * 0.5); run();

        double avg = getExpectedAverageDegree();

        if (avg > avgDeg) {
            average_r = avg;
            gamma_r = gamma;
        } else {
            average_l = avg;
            gamma_l = gamma;
        }
    }

    if (avgDeg - average_l < average_r - avgDeg) {
        setGamma(gamma_l);
    } else {
        setGamma(gamma_r);
    }
}


NetworKit::count NetworKit::PowerlawDegreeSequence::getMinimumDegree() const {
    return minDeg;
}

void NetworKit::PowerlawDegreeSequence::run() {
    cumulativeProbability.clear();
    cumulativeProbability.reserve(maxDeg - minDeg + 1);

    double sum = 0;

    for (count d = maxDeg; d >= minDeg; --d) {
        sum += std::pow(static_cast<double>(d), gamma);
        cumulativeProbability.push_back(sum);
    }

    for (double & prob : cumulativeProbability) {
        prob /= sum;
    }

    cumulativeProbability.back() = 1.0;

    hasRun= true;
}

double NetworKit::PowerlawDegreeSequence::getExpectedAverageDegree() const {
    assureFinished();
    double average = cumulativeProbability[0] * maxDeg;
    for (count i = 1; i < cumulativeProbability.size(); ++i) {
        average += (cumulativeProbability[i] - cumulativeProbability[i-1]) * static_cast<double>(maxDeg - i);
    }

    return average;
}

std::vector< NetworKit::count > NetworKit::PowerlawDegreeSequence::getDegreeSequence(NetworKit::count numNodes) const {
    std::vector<count> degreeSequence;

    assureFinished();
    degreeSequence.reserve(numNodes);
    count degreeSum = 0;

    for (count i = 0; i < numNodes; ++i) {
        degreeSequence.push_back(getDegree());
        degreeSum += degreeSequence.back();
    }

    if (degreeSum % 2 != 0) {
        (*std::max_element(degreeSequence.begin(), degreeSequence.end()))--;
    }

    return degreeSequence;
}

NetworKit::count NetworKit::PowerlawDegreeSequence::getDegree() const {
    assureFinished();
    return maxDeg - std::distance(cumulativeProbability.begin(), std::lower_bound(cumulativeProbability.begin(), cumulativeProbability.end(), Aux::Random::probability()));
}
