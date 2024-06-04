// -*- LSST-C++ -*-
/*
 * This file is part of gauss2d.
 *
 * Developed for the LSST Data Management System.
 * This product includes software developed by the LSST Project
 * (https://www.lsst.org).
 * See the COPYRIGHT file at the top-level directory of this distribution
 * for details of code ownership.
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <https://www.gnu.org/licenses/>.
 */

#ifndef LSST_GAUSS2D_TYPE_NAME_H
#define LSST_GAUSS2D_TYPE_NAME_H

#include <string>
#include <string_view>

namespace lsst::gauss2d {

template <typename T>
constexpr std::string_view type_name();

template <>
constexpr std::string_view type_name<void>() {
    return "void";
}

namespace detail {

// Adapted from https://stackoverflow.com/a/64490578

using type_name_prober = void;

template <typename T>
constexpr std::string_view wrapped_type_name() {
#ifdef __clang__
    return __PRETTY_FUNCTION__;
#elif defined(__GNUC__)
    return __PRETTY_FUNCTION__;
#elif defined(_MSC_VER)
    return __FUNCSIG__;
#else
#error "Unsupported compiler"
#endif
}

constexpr std::size_t wrapped_type_name_prefix_length() {
    return wrapped_type_name<type_name_prober>().find(type_name<type_name_prober>());
}

constexpr std::size_t wrapped_type_name_suffix_length() {
    return wrapped_type_name<type_name_prober>().length() - wrapped_type_name_prefix_length()
           - type_name<type_name_prober>().length();
}

constexpr std::string_view NAMESPACE_SEPARATOR = "::";
constexpr auto NAMESPACE_SEPARATOR_LEN = NAMESPACE_SEPARATOR.size();

}  // namespace detail

/**
 * Get a string representation of an arbitrary C++ type.
 *
 * @tparam T The type to stringify
 * @return A string representation of the type's name
 *
 * @note Adapted from https://stackoverflow.com/a/64490578
 */
template <typename T>
constexpr std::string_view type_name() {
    constexpr auto wrapped_name = detail::wrapped_type_name<T>();
    constexpr auto prefix_length = detail::wrapped_type_name_prefix_length();
    constexpr auto suffix_length = detail::wrapped_type_name_suffix_length();
    constexpr auto type_name_length = wrapped_name.length() - prefix_length - suffix_length;
    return wrapped_name.substr(prefix_length, type_name_length);
}

/**
 * Get a string representation of an arbitrary C++ type, potentially modifying
 * its namespace prefix.
 *
 * @tparam T The type to stringify.
 * @param strip_namespace Whether to strip the namespace prefix entirely.
 * @param namespace_str A string to replace the standard C++ namespace
 *      separator (i.e. ::) with; generally . for Python.
 * @return A string representation of the type's name, with modified namespace
 *      prefix.
 */
template <typename T>
std::string type_name_str(bool strip_namespace = false,
                          std::string_view namespace_str = detail::NAMESPACE_SEPARATOR) {
    std::string name = std::string(type_name<T>());
    if (strip_namespace) {
        return name.substr(name.find_last_of(':') + 1, std::string::npos);
    } else if (namespace_str != detail::NAMESPACE_SEPARATOR) {
        auto pos = name.find(detail::NAMESPACE_SEPARATOR, 0);
        const auto n_replace = namespace_str.size();
        while (pos != std::string::npos) {
            name.replace(pos, detail::NAMESPACE_SEPARATOR_LEN, namespace_str);
            pos += n_replace;
            pos = name.find(detail::NAMESPACE_SEPARATOR, pos);
        }
    }
    return name;
}

}  // namespace lsst::gauss2d
#endif
