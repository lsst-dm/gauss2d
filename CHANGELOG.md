# Change Log

All notable changes to this project will be documented in this file.
This project adheres to [Semantic Versioning](https://semver.org/).

### [0.1.1] 2024-07-19

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

[0.1.1]: https://github.com/lsst-dm/gauss2d/compare/0.1.0...0.1.1
[0.1.0]: https://github.com/lsst-dm/gauss2d/compare/53bc2990d...0.1.0
