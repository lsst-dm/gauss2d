#ifndef LSST_GAUSS2D_TO_STRING_H
#define LSST_GAUSS2D_TO_STRING_H

#include <algorithm>
#include <memory>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace lsst::gauss2d {

template <typename T>
std::string to_string_float(const T value, const int precision = 6, const bool scientific = true) {
    std::ostringstream out;
    out.precision(precision);
    if (scientific) {
        out << std::scientific;
    } else {
        out << std::fixed;
    }
    out << value;
    return std::move(out).str();
}

template <template <typename...> class Container, class Value>
std::string to_string_float_iter(const Container<Value>& container, const int precision = 6,
                                 const bool scientific = true) {
    std::string str = "[";
    for (const auto& value : container) {
        str += to_string_float(value, precision, scientific) + ", ";
    }
    return str.substr(0, str.size() - 2 * (container.size() > 0)) + "]";
}

template <template <typename...> class Container, class Value>
std::string to_string_iter(const Container<Value>& container) {
    std::string str = "[";
    for (const auto& value : container) {
        str += std::to_string(value) + ", ";
    }
    return str.substr(0, str.size() - 2 * (container.size() > 0)) + "]";
}

}  // namespace lsst::gauss2d

#endif
