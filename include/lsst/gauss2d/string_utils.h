#ifndef LSST_GAUSS2D_STRING_UTILS_H
#define LSST_GAUSS2D_STRING_UTILS_H

#include <algorithm>
#include <memory>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace lsst::gauss2d {

/**
 * Replace a token inside a target string with another string.
 *
 * @param target The string to replace and return.
 * @param token The token to replace.
 * @param replacement The token to replace with.
 * @return The target with replacements, if any.
 */
template <typename T, typename R>
std::string replace_all(std::string target, T token, R replacement) {
    auto pos = target.find(token, 0);
    const auto n_token = token.size();
    const auto n_replace = replacement.size();
    while (pos != std::string::npos) {
        target.replace(pos, n_token, replacement);
        pos += n_replace;
        pos = target.find(token, pos);
    }
    return target;
}

/**
 * Replace a token inside a target string with nothing.
 *
 * @param target The string to replace and return.
 * @param token The token to replace with an empty string.
 * @return The target with replacements, if any.
 */
template <typename T>
std::string replace_all_none(std::string target, T token) {
    auto pos = target.find(token, 0);
    const auto n_token = token.size();
    while (pos != std::string::npos) {
        target.replace(pos, n_token, "");
        pos = target.find(token, pos);
    }
    return target;
}

template std::string replace_all<const std::string, const std::string>(std::string target,
                                                                       const std::string token,
                                                                       const std::string replacement);
template std::string replace_all<std::string_view, const std::string>(std::string target,
                                                                      std::string_view token,
                                                                      const std::string replacement);
template std::string replace_all<const std::string, std::string_view>(std::string target,
                                                                      const std::string token,
                                                                      std::string_view replacement);
template std::string replace_all<std::string_view, std::string_view>(std::string target,
                                                                     std::string_view token,
                                                                     std::string_view replacement);

template std::string replace_all_none<const std::string>(std::string target, const std::string token);
template std::string replace_all_none<std::string_view>(std::string target, std::string_view token);
}  // namespace lsst::gauss2d

#endif
