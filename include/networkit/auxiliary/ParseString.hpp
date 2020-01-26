/*
 * ParseString.hpp
 *
 * Created: 2019-07-20
 * Author: Armin Wiebigke
 */

#ifndef NETWORKIT_AUXILIARY_PARSE_STRING_HPP_
#define NETWORKIT_AUXILIARY_PARSE_STRING_HPP_

#include <string>

namespace Aux {

/**
 * Convert a string to a floating point number.
 * @param str
 * @param decimalChar The character representing the decimal point, either ',' or '.'
 * @return
 */
double stringToDouble(std::string str, char decimalChar = '.');

}

#endif // NETWORKIT_AUXILIARY_PARSE_STRING_HPP_
