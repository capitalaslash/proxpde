# ProXPDE - Prototyping with eXpression Templates for Partial Differential Equations

Templated c++-20 implementation (mainly) of the finite element method for PDE solving.

## Installation

The recommended way to install all the dependencies is via `spack`.

Install `spack` following the instruction from
[readthedocs](https://spack.readthedocs.io/en/latest/). Use the spack environment
provided in the main source directory

```
spack env activate spack_env
spack install
```

to install all the required dependencies. After that, just activating the environment
will set the terminal for compilation and execution.

### Dependencies

* `Eigen` (minimum version 3.4)

* `libfmt`

* `hdf5`

* `OpenMP`

* `pugixml`

* `UMFPack`

* `yaml-cpp`
