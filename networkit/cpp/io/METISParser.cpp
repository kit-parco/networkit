// no-networkit-format
/*
 * METISParser.cpp
 *
 *  Created on: 27.11.2012
 *      Author: Christian Staudt
 */

#include <networkit/auxiliary/NumberParsing.hpp>
#include <networkit/auxiliary/Enforce.hpp>
#include <networkit/auxiliary/Log.hpp>
#include <networkit/io/METISParser.hpp>

#include <stdexcept>

namespace NetworKit {

/**
 * Extract a vector of indices from a line in the file.
 *
 * @param[in]  line    line from input file containing node indices
 * @param[in]  ignoreFirst  number of values to ignore [in case the METIS file contains node weightes]
 *
 * @param[out]  adjacencies  node indices extracted from line
 */
static std::vector<node> parseLine(const std::string& line, count ignoreFirst = 0) {
    auto it = line.begin();
    auto end = line.end();
    std::vector<node> adjacencies;
    node v;
    index i = 0;
    DEBUG(ignoreFirst);
    while (i < ignoreFirst) {
        // parse first values but ignore them.
        double dummy;
        std::tie(dummy, it) = Aux::Parsing::strTo<double>(it,end);
        DEBUG("ignored: ",dummy);
        ++i;
    }
    while(it != end) {
        std::tie(v, it) = Aux::Parsing::strTo<node>(it,end);
        adjacencies.push_back(v);
    }
    return adjacencies;
}

/**
 * Extract a vector of indices and edge weights from a line in the file.
 *
 * @param[in]  line    line from input file containing node indices
 * @param[in]  ignoreFirst  number of values to ignore [in case the METIS file contains node weightes]
 *
 * @param[out]  adjacencies  node indices including edge weights extracted from line
 */
static std::vector<std::pair<node,double>> parseWeightedLine(std::string line, count ignoreFirst = 0) {
    auto it = line.begin();
    auto end = line.end();
    std::vector<std::pair<node,double>> adjacencies;
    node v;
    double weight;
    index i = 0;
    DEBUG(ignoreFirst);
    std::stringstream content;
    while (i < ignoreFirst) {
        // parse first values but ignore them.
        double dummy;
        std::tie(dummy, it) = Aux::Parsing::strTo<double>(it,end);
        DEBUG("ignored: ",dummy);
        ++i;
    }
    while(it != end) {
        try {
            std::tie(v, it) = Aux::Parsing::strTo<node>(it,end);
            std::tie(weight, it) = Aux::Parsing::strTo<double,decltype(it),Aux::Checkers::Enforcer>(it,end);
            adjacencies.emplace_back(v, weight);
            content << v << " " << weight << "\t";
        } catch (const std::exception &e) {
            ERROR("malformed line; not all edges have been read correctly");
            break;
        }
    }
    DEBUG(content.str());
    return adjacencies;
}

METISParser::METISParser(const std::string &path) : graphFile(path) {
    if (!(this->graphFile)) {
        ERROR("invalid graph file: " , path);
        throw std::runtime_error("invalid graph file");
    }
}

std::tuple<count, count, index, count> METISParser::getHeader() {
    // handle header line
    count n;    // number of nodes
    count m;    // number of edges
    index fmt = 0;    // weighted or unweighted graph
    count ncon = 0;    // number of node weights

    std::string line = "";
    Aux::enforceOpened(this->graphFile);

    if (std::getline(this->graphFile, line)) {
        // ignore comment lines
        while (line[0] == '%') {
            std::getline(this->graphFile, line);
        }

        std::vector<uint64_t> tokens = parseLine(line);
        n = tokens[0];
        m = tokens[1];
        if (tokens.size() == 2) {
            return std::tuple<count, count, index, count>(n,m,fmt,ncon);
        }
        if (tokens.size() >= 3) {
            fmt = tokens[2];
            if (fmt >= 2) {
                WARN("nodes are weighted; node weights will be ignored");
            }
            if (tokens.size() == 4) {
                ncon = tokens[3];
            } else {
                ncon = 1;
            }
        }
        return std::tuple<count, count, index, count>(n,m,fmt,ncon);
    } else {
        ERROR("getline not successful");
        throw std::runtime_error("getting METIS file header failed");
        return std::tuple<count, count, index, count>(0,0,0,0);
    }
}

bool METISParser::hasNext() {
    // if graph file has lines left, return true
    return this->graphFile.good();
}

std::vector<node> METISParser::getNext(count ignoreFirst) {

    std::string line;
    do {
        std::getline(this->graphFile, line);
        // check for comment line starting with '%'
        if (line[0] == '%') {
            continue;
        } else {
            return parseLine(line,ignoreFirst);
        }
    } while (true);

    throw std::runtime_error("bad METIS file structure");
}

std::vector<std::pair<node,double>> METISParser::getNextWithWeights(count ignoreFirst) {

    std::string line;
    do {
        std::getline(this->graphFile, line);
        // check for comment line starting with '%'
        if (line[0] == '%') {
            continue;
        } else {
            return parseWeightedLine(line,ignoreFirst);
        }
    } while (true);

    throw std::runtime_error("bad METIS file structure");
}

} /* namespace NetworKit */
