// no-networkit-format
/*
 *  MatrixMarketReader.cpp
 *
 *  Created on: 25.07.2014
 *      Author: dhoske
 */

#include <fstream>
#include <stdexcept>
#include <string>

#include <networkit/io/MatrixMarketReader.hpp>

namespace NetworKit {

namespace {
    static constexpr char COMMENT_CHAR = '%';
    static const std::string MAGIC = "%%matrixmarket";

    std::string tolower(const std::string& str) {
        std::string out;
        std::transform(str.begin(), str.end(), std::back_inserter(out), ::tolower);
        return out;
    }
}

CSRMatrix MatrixMarketReader::read(const std::string &path) {
    std::ifstream in(path);
    if (!in.is_open()) {
        throw std::runtime_error("could not open: " + path);
    }
    return read(in);
}

CSRMatrix MatrixMarketReader::read(std::istream& in) {
    enum { FIRST_LINE, HEADER, ENTRIES } state = FIRST_LINE;

    count nrows = 0, ncols = 0, nzeroes;
    std::vector<Triplet> triplets;
    std::string line;
    bool weighted = true;
    bool symmetric = false;
    while (getline(in, line)) {
        std::istringstream in_line(line);

        if (state == FIRST_LINE) {
            // First line: magic <keyword>*

            // Read id of file format at start of file
            std::string magic, object, format, data, qualifier;
            in_line >> magic >> object >> format >> data >> qualifier;
            magic = tolower(magic);
            object = tolower(object);
            format = tolower(format);
            data = tolower(data);
            qualifier = tolower(qualifier);
            if (magic != MAGIC) {
                throw std::runtime_error(MAGIC + " not found");
            }
            if (format != "coordinate") {
                throw std::runtime_error("Unsupported format: " + format);
            }
            if (data == "pattern") {
                weighted = false;
            } else if (data != "real") {
                throw std::runtime_error("Unsupported data type: " + data);
            }

            if (qualifier == "symmetric") {
                symmetric = true;
            }

            state = HEADER;
        } else if (line.empty() || line[0] == COMMENT_CHAR) {
            // Skip line if empty or comment
            continue;
        } else if (state == HEADER) {
            // Header: nrows, ncols, nzeroes
            in_line >> nrows >> ncols >> nzeroes;
            if (in_line.fail()) {
                throw std::runtime_error("expected three non-negative integers in header line:" + line);
            }

            state = ENTRIES;
        } else {
            // Entry: row col value
            count i, j;
            double val;
            in_line >> i >> j;
            i--;
            j--;
            if (weighted) {
                in_line >> val;
            } else {
                val = 1.0;
            }

            if (in_line.fail()) {
                throw std::runtime_error("expected three entries in entry line: " + line);
            }
            if (i >= nrows) {
                throw std::runtime_error("invalid index: " + std::to_string(i));
            }
            if (j >= ncols) {
                throw std::runtime_error("invalid index: " + std::to_string(j));
            }

            triplets.push_back({i,j,val});

            if (i != j && symmetric) {
                triplets.push_back({j,i,val});
            }
        }

    }

    return CSRMatrix(nrows, ncols, triplets);
}

} // namespace NetworKit
