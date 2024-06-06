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

#ifndef LSST_GAUSS2D_OBJECT_H
#define LSST_GAUSS2D_OBJECT_H

#include <string>

namespace lsst::gauss2d {

/**
 * A generic object from the gauss2d library.
 *
 * Objects have string representations that are guaranteed to be valid in C++
 * and should also be valid in Python bindings, if implemented correctly.
 * The interface may be expanded in the future.
 */
class Object {
public:
    /// The C++ namespace separator
    static constexpr std::string_view CC_NAMESPACE_SEPARATOR = "::";
    static constexpr std::string_view NULL_STR_GENERAL = "None";
    static constexpr std::string_view PY_NAMESPACE_SEPARATOR = ".";

    static std::string_view null_str(const std::string_view &namespace_separator) {
        return namespace_separator == CC_NAMESPACE_SEPARATOR ? "nullptr" : NULL_STR_GENERAL;
    }

    /**
     * Return a full, callable string representation of this.
     *
     * @param name_keywords Whether to prefix arguments with "{name}=",
     *      where name is the arg name in the header (as with keyword
     *      arguments in Python).
     * @param namespace_separator The string to use to delimit namespaces,
     *      i.e. :: in C++ and . in Python.
     * @return
     *      A callable string representation of this, which should return an
     *      an identical object to this.
     *
     * @note The representation with name_keywords=false must be callable
     * in C++. The representation with name_keywords=true should be callable
     * in Python, if there are any bindings.
     */
    virtual std::string repr(bool name_keywords = false,
                             std::string_view namespace_separator = CC_NAMESPACE_SEPARATOR) const
            = 0;
    /// Return a brief, human-readable string representation of this.
    virtual std::string str() const = 0;

    friend std::ostream &operator<<(std::ostream &out, const Object &obj) {
        out << obj.str();
        return out;
    }

    virtual ~Object() = default;
};

template <typename T>
std::string repr_ptr(T ptr, bool name_keywords, std::string_view namespace_separator) {
    return ptr ? ptr->repr(name_keywords, namespace_separator)
               : std::string(Object::null_str(namespace_separator));
}

template <typename T>
std::string repr_iter_ptr(const T &container, bool name_keywords = false,
                          std::string_view namespace_separator = Object::CC_NAMESPACE_SEPARATOR) {
    std::string str = "[";
    for (auto &obj : container) {
        str += repr_ptr(obj, name_keywords, namespace_separator) + ", ";
    }
    auto size_str = str.size();
    if (size_str > 1) str = str.substr(0, size_str - 2);
    return str + "]";
}

template <typename T, bool is_wrapper>
std::string repr_iter_ref(const T &container, bool name_keywords = false,
                          std::string_view namespace_separator = Object::CC_NAMESPACE_SEPARATOR) {
    std::string str = "[";
    for (const auto &obj : container) {
        str += (is_wrapper ? obj.get() : obj).repr(name_keywords, namespace_separator) + ", ";
    }
    auto size_str = str.size();
    if (size_str > 1) str = str.substr(0, size_str - 2);
    return str + "]";
}

template <typename T>
std::string str_ptr(T ptr) {
    return ptr ? ptr->str() : std::string(Object::NULL_STR_GENERAL);
}

template <typename T>
std::string str_iter_ptr(const T &container) {
    std::string str = "[";
    for (const auto &obj : container) {
        str += str_ptr(obj) + ", ";
    }
    auto size_str = str.size();
    if (size_str > 1) str = str.substr(0, size_str - 2);
    return str + "]";
}

template <typename T, bool is_wrapper>
std::string str_iter_ref(const T &container) {
    std::string str = "[";
    for (const auto &obj : container) {
        str += (is_wrapper ? obj.get() : obj).str() + ", ";
    }
    auto size_str = str.size();
    if (size_str > 1) str = str.substr(0, size_str - 2);
    return str + "]";
}

}  // namespace lsst::gauss2d
#endif
