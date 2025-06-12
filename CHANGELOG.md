# Change Log

All notable changes to this project will be documented in this file.
This project adheres to [Semantic Versioning](https://semver.org/).

### [0.2.1] 2025-06-10

* Fixed: Removed a redundant build-cc call in EUPS config
* See [DM-45703](https://rubinobs.atlassian.net/browse/DM-45703) for details.

### [0.2.0] 2025-04-29

* Changed: Updated default C++ standard to C++20 from C++17
* See [DM-50434](https://rubinobs.atlassian.net/browse/DM-50434) for details.

### [0.1.3] 2024-07-30

* Changed: Updated Python module README build instructions
* Added: Instructions on adding .pth file for standalone install in README
* Added: PyImage GaussianEvaluator test case
* Fixed: segfault in pybind11 module on Pybind 3.12
* Fixed: Made evaluator in make_gaussians_pixel a unique_ptr, not shared
* Fixed: Added missing override
* See [DM-45473](https://rubinobs.atlassian.net/browse/DM-45473) for details.

### [0.1.2] 2024-07-23

* Fixed: Added pytest.ini_options to python/pyproject.toml
* Fixed: Added missing LSST_LIBRARY_PATH entry in ups/gauss2d.table
* This is a hotfix for [DM-45346](https://rubinobs.atlassian.net/browse/DM-45346).

### [0.1.1] 2024-07-22

* Added: doctest upgrade to 2.4.11
* Added: clean() command in ups/eupspkg.cfg.sh.
* Fixed: Added missing override in EllipseValues::set_h
* Fixed: Only test Image str/repr on linux (to be investigated in [DM-45308](https://rubinobs.atlassian.net/browse/DM-45308))
* Fixed: lib64 refs in python/meson.build and meson.build replaced with lib.
* Fixed: --libdir "lib" added to eups meson config command.
* See [DM-45346](https://rubinobs.atlassian.net/browse/DM-45346) for details.

### [0.1.0] 2024-07-09

* Changed: Initial release, forked to https://github.com/lsst/gauss2d.
* See [DM-43906](https://rubinobs.atlassian.net/browse/DM-43906) for details. 

[0.2.1]: https://github.com/lsst-dm/gauss2d/compare/0.2.0...0.2.1
[0.2.0]: https://github.com/lsst-dm/gauss2d/compare/0.1.3...0.2.0
[0.1.3]: https://github.com/lsst-dm/gauss2d/compare/0.1.2...0.1.3
[0.1.2]: https://github.com/lsst-dm/gauss2d/compare/0.1.1...0.1.2
[0.1.1]: https://github.com/lsst-dm/gauss2d/compare/0.1.0...0.1.1
[0.1.0]: https://github.com/lsst-dm/gauss2d/compare/53bc2990d...0.1.0
