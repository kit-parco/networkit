// no-networkit-format
/*
 * DGSStreamParser.cpp
 *
 *  Created on: 23.12.2013
 *      Author: cls
 */

#include <networkit/auxiliary/Log.hpp>
#include <networkit/auxiliary/StringTools.hpp>
#include <networkit/dynamics/DGSStreamParser.hpp>

namespace NetworKit {

DGSStreamParser::DGSStreamParser(const std::string &path, bool mapped, node baseIndex)
    : dgsFile(path), mapped(mapped), baseIndex(baseIndex), nextNode(0) {}

std::vector<GraphEvent> DGSStreamParser::getStream() {
    if (!dgsFile.is_open()) {
        throw std::runtime_error("DGS input file could not be opened.");
    }

    std::vector<GraphEvent> stream; // stream containing the events

    std::string line;
    count lineCount = 0;

    std::string cookie = "DGS004";
    std::getline(dgsFile, line); // get DGS version
    ++lineCount;
    if (line.compare(0, cookie.size(), cookie)) {
        ERROR("found " , line , " instead of " , cookie , " in first line");
        throw std::runtime_error("expected cookie in first line");
    }

    std::getline(dgsFile, line);
    ++lineCount;
    INFO("DGS stream description: ", line);

    if (mapped) {
        // mapped format

        /**
         * Maps key string to consecutive, 0-based node id.
         */
        auto map = [&](const std::string &key) {
            auto iter = this->key2id.find(key);
            if (iter == key2id.end()) {
                key2id[key] = nextNode;
                nextNode++;
                return key2id[key];
            } else {
                return iter->second;
            }
        };

        while (std::getline(dgsFile, line)) {
            ++lineCount;
            std::vector<std::string> split = Aux::StringTools::split(line);
            std::string tag = split[0];

            // parse commands
            if (tag.compare("st") == 0) { // clock
                stream.emplace_back(GraphEvent::TIME_STEP);
            } else if (tag.compare("an") == 0) { // add node
                node u = map(split[1]);
                auto ev = GraphEvent(GraphEvent::NODE_ADDITION, u);
                stream.push_back(ev);
            } else if (tag.compare("ae") == 0) { // add edge
                node u = map(split[2]);
                node v = map(split[3]);
                edgeweight w = 1.0;
                if (split.size() >= 5) {
                    w = std::stod(Aux::StringTools::split(split[4], '=')[1]); // weight=<w>
                }
                auto ev = GraphEvent(GraphEvent::EDGE_ADDITION, u, v, w);
                stream.push_back(ev);
            } else if (tag.compare("ce") == 0) {
                // update edge. Only the "weight" attribute is supported so far
                std::vector<std::string> uvs = Aux::StringTools::split(split[1], '-');
                node u = map(uvs[0]);
                node v = map(uvs[1]);
                edgeweight w = std::stod(Aux::StringTools::split(split[2], '=')[1]); // weight=<w>
                auto ev = GraphEvent(GraphEvent::EDGE_WEIGHT_UPDATE, u, v, w);
                stream.push_back(ev);
            } else if (tag.compare("ie") == 0) { // update edge. Only the "weight" attribute is supported so far
                std::vector<std::string> uvs = Aux::StringTools::split(split[1], '-');
                node u = map(uvs[0]);
                node v = map(uvs[1]);
                edgeweight w = std::stod(Aux::StringTools::split(split[2], '=')[1]); // weight=<w>
                auto ev = GraphEvent(GraphEvent::EDGE_WEIGHT_INCREMENT, u, v, w);
                stream.push_back(ev);
            } else if (tag.compare("de") == 0) {
                std::vector<std::string> uvs = Aux::StringTools::split(split[1], '-');
                node u = map(uvs[0]);
                node v = map(uvs[1]);
                auto ev = GraphEvent(GraphEvent::EDGE_REMOVAL, u, v);
                stream.push_back(ev);
            } else if (tag.compare("dn") == 0) {
                node u = map(split[1]);
                auto ev = GraphEvent(GraphEvent::NODE_REMOVAL, u);
                stream.push_back(ev);
            } else if (tag.compare("rn") == 0) {
                node u = map(split[1]);
                auto ev = GraphEvent(GraphEvent::NODE_RESTORATION, u);
                stream.push_back(ev);
            } else {
                ERROR("malformed line (" , lineCount , ") : " , line);
                throw std::runtime_error("malformed line in .DGS file");
            }
        }

    } else {
        // direct format
        auto offset = [&](node u) {
            return (u - baseIndex);
        };

        while (std::getline(dgsFile, line)) {
            ++lineCount;
            std::vector<std::string> split = Aux::StringTools::split(line);
            std::string tag = split[0];

            // parse commands
            if (tag.compare("st") == 0) { // clock
                stream.emplace_back(GraphEvent::TIME_STEP);
            } else if (tag.compare("an") == 0) { // add node
                node u = offset(std::stoul(split[1]));
                stream.emplace_back(GraphEvent::NODE_ADDITION, u);
            } else if (tag.compare("ae") == 0) { // add edge
                node u = offset(std::stoul(split[2]));
                node v = offset(std::stoul(split[3]));
                edgeweight w = std::stod(Aux::StringTools::split(split[4], '=')[1]); // weight=<w>
                stream.emplace_back(GraphEvent::EDGE_ADDITION, u, v, w);
            } else if (tag.compare("ce") == 0) {
                // update edge. Only the "weight" attribute is supported so far
                std::vector<std::string> uvs = Aux::StringTools::split(split[1], '-');
                node u = offset(std::stoul(uvs[0]));
                node v = offset(std::stoul(uvs[1]));
                edgeweight w = std::stod(Aux::StringTools::split(split[2], '=')[1]); // weight=<w>
                stream.emplace_back(GraphEvent::EDGE_WEIGHT_UPDATE, u, v, w);
            } else if (tag.compare("ie") == 0) {
                // update edge. Only the "weight" attribute is supported so far
                std::vector<std::string> uvs = Aux::StringTools::split(split[1], '-');
                node u = offset(std::stoul(uvs[0]));
                node v = offset(std::stoul(uvs[1]));
                edgeweight w = std::stod(Aux::StringTools::split(split[2], '=')[1]); // weight=<w>
                stream.emplace_back(GraphEvent::EDGE_WEIGHT_INCREMENT, u, v, w);
            } else if (tag.compare("de") == 0) {
                std::vector<std::string> uvs = Aux::StringTools::split(split[1], '-');
                node u = offset(std::stoul(uvs[0]));
                node v = offset(std::stoul(uvs[1]));
                stream.emplace_back(GraphEvent::EDGE_REMOVAL, u, v);
            } else if (tag.compare("dn") == 0) {
                node u = offset(std::stoul(split[1]));
                stream.emplace_back(GraphEvent::NODE_REMOVAL, u);
            } else if (tag.compare("rn") == 0) {
                node u = offset(std::stoul(split[1]));
                stream.emplace_back(GraphEvent::NODE_RESTORATION, u);
            } else {
                ERROR("malformed line (" , lineCount , ") : " , line);
                throw std::runtime_error("malformed line in .DGS file");
            }
        }
    }

    return stream;
}

} /* namespace NetworKit */
