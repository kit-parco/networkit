/*
 * ParseString.cpp
 *
 * Created: 2019-07-20
 * Author: Armin Wiebigke
 */

#include <algorithm>
#include <stdexcept>
#include <cassert>

#include <networkit/auxiliary/ParseString.hpp>

namespace Aux {

void swapChars(std::string &str, char char1, char char2, char tmpChar);

double stringToDouble(std::string str, char decimalChar) {
    assert(decimalChar == '.' || decimalChar == ',');
    // Use decimal separator of current locale
    if ((decimalChar == '.' && std::stod("0.3") == 0) ||
        (decimalChar == ',' && std::stod("0,3") == 0)) {
        swapChars(str, '.', ',', '#');
    }
    try {
        return std::stod(str);
    } catch (const std::invalid_argument &e) {
        throw std::invalid_argument("Could not convert the string \"" + str + "\" into a double.");
    }

}

void swapChars(std::string &str, char char1, char char2, char tmpChar) {
    assert(str.find(tmpChar) == std::string::npos);
    std::replace(str.begin(), str.end(), char1, tmpChar);
    std::replace(str.begin(), str.end(), char2, char1);
    std::replace(str.begin(), str.end(), tmpChar, char2);
}

}