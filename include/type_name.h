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

#ifndef GAUSS2D_TYPE_NAME_H
#define GAUSS2D_TYPE_NAME_H

#include <string_view>

namespace gauss2d {

// https://stackoverflow.com/questions/81870/is-it-possible-to-print-a-variables-type-in-standard-c/64490578#64490578
template <typename T>
constexpr std::string_view type_name();

template <>
constexpr std::string_view type_name<void>() {
    return "void";
}

namespace detail {

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

}  // namespace detail

template <typename T>
constexpr std::string_view type_name() {
    constexpr auto wrapped_name = detail::wrapped_type_name<T>();
    constexpr auto prefix_length = detail::wrapped_type_name_prefix_length();
    constexpr auto suffix_length = detail::wrapped_type_name_suffix_length();
    constexpr auto type_name_length = wrapped_name.length() - prefix_length - suffix_length;
    return wrapped_name.substr(prefix_length, type_name_length);
}

}  // namespace gauss2d
#endif