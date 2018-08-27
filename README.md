**This package is deprecated and will not be maintained anymore because all use cases are covered by the [DifferentialEquations.jl](http://docs.juliadiffeq.org/latest/index.html) ecosystem.** 

# Dopri.jl

[![Build Status](https://travis-ci.org/helgee/Dopri.jl.svg)](https://travis-ci.org/helgee/Dopri.jl)
[![Build status](https://ci.appveyor.com/api/projects/status/aj34gxo8noq5lecp?svg=true)](https://ci.appveyor.com/project/helgee/dopri-jl)
[![Coverage Status](https://coveralls.io/repos/helgee/Dopri.jl/badge.svg?branch=master&service=github)](https://coveralls.io/github/helgee/Dopri.jl?branch=master)

Dopri.jl is a Julia wrapper for the DOPRI5 and DOP853 integrators by [Ernst Hairer](http://www.unige.ch/~hairer/software.html).

## Install

Dopri.jl can then be installed through Julia's package manager.
On Linux and macOS you need to have either __Gfortran__ or the __Intel Fortran Compiler__ installed to be able to build the binary dependencies.

```julia
Pkg.add("Dopri")
```

## Usage

The API is modelled after that of [ODE.jl](https://github.com/JuliaLang/ODE.jl).

```julia
tout, yout = dopri5(F!, y0, tspan; keywords...)
tout, yout = dop853(F!, y0, tspan; keywords...)
```

The only differences are the mutating user function `F!` (see below) and a few additional keyword arguments.

The following keyword arguments are supported:

* `atol`: Vector of absolute tolerances. Default: `sqrt(eps())`.
* `rtol`: Vector of realtive tolerances. Default: `1e-6`.
* `points`
    * `= :all`: Output is given for each value in `tspan` and all intermediate solver steps.
    * `= :specified`: Output is given for each value in `tspan`.
    * `= :last`: Output is given for the final step only.
* `solout`: User function that is called after every successfull integration step (see below).
* `dense`: Vector of indices for which dense output shall be performed. Default: Dense output is provided for all components. Is also set automatically for `points=:specified` and `points=:all`.
* `verbose`: Print messages from `DOPRI5` and `DOP853`. Default: `false`.

### User Functions
```julia
F!(f, t, y) = ...
```

* `f`: `f=dy/dt`
* `t`: Current step.
* `y`: Current state vector.

```julia
solout!(told, t, y, contd) = ... return dopricode[:nominal]
```

If the user supplies a `solout` function, it will be called after every successful integration step. Within `solout` the `contd` function can be used to approximate state vector components between the current and the preceding integration step via dense output.

* `told`: Last step.
* `t`: Current step.
* `y`: Current state vector.
* `contd`: `yi = contd(i, t1)` Get dense output `yi` for component `i` at `t1` with `told < t1 < t`.

`solout` must return one of the following return codes:

* `dopricode[:nominal]`: If the integration shall continue nominally.
* `dopricode[:altered]`: If numerical solution was altered in `solout`.
* `dopricode[:abort]`: If the integration shall be stopped.

### Passing Additional Parameters to User Functions

A closure can be used to pass additional parameters to the user functions, e.g.:

```julia
parameter = pi
tout, yout = dop853((f, t, y) -> F!(f, t, y, parameter), y0, tspan)
```
