# This file is part of gauss2d.
#
# Developed for the LSST Data Management System.
# This product includes software developed by the LSST Project
# (https://www.lsst.org).
# See the COPYRIGHT file at the top-level directory of this distribution
# for details of code ownership.
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <https://www.gnu.org/licenses/>.

import lsst.gauss2d as g2d


def test_bindings():
    errors = []
    bad_args_msg = "incompatible function arguments"

    for name_class in dir(g2d):
        klass = getattr(g2d, name_class)

        if callable(klass):
            try:
                obj = klass()
            except ValueError:
                # Likely non-trivial __init__
                continue
            except TypeError as e:
                msg = str(e)
                if not (
                    msg.endswith("No constructor defined!")
                    or msg.startswith("__init__(): incompatible constructor arguments.")
                    or msg.startswith("extend_path() missing")
                    or (msg.startswith("make_gaussians_pixel") and (bad_args_msg in msg))
                ):
                    errors.append(e)
                    # str and repr should not raise
            try:
                repr(obj)
            except Exception as e:
                errors.append(e)
            try:
                str(obj)
            except Exception as e:
                errors.append(e)

    if errors:
        print('\n'.join((str(e) for e in errors)))
        assert not 'errors not empty; see stdout'
