#  JuMP Demo at Dagstuhl Seminar 18081

This site contains materials for a demonstration of [JuMP](https://github.com/JuliaOpt/JuMP.jl) at the [Dagstuhl Seminar 18081: Designing and Implementing Algorithms for Mixed-Integer Nonlinear Optimization](http://www.dagstuhl.de/en/program/calendar/semhp/?semnr=18081)

## Install Julia

You should use the latest version of Julia v0.6.2.
Binaries of Julia for all platforms are available [here](http://julialang.org/downloads/).


## Install Basic JuMP packages

To start off, we will be using the following packages:
- JuMP
- Clp/Cbc

To install them, run the following code:
```julia
Pkg.add("JuMP")
Pkg.add("Clp")
```
If you have a previous installation of Julia,
be sure to update your packages to the latest version by running ``Pkg.update()``.

To test that your installation is working, run the following code (the first time you run the code you may see the message "INFO: Precompiling module JuMP" for a few minutes):

```julia
using JuMP
m = Model()
@variable(m, x >= 0)
@variable(m, y >= 0)
@constraint(m, 2x + y <= 1)
@objective(m, Max, x+y)
status = solve(m)
@show status
@show getvalue(y)
```

The output should be:

```
status = :Optimal
getvalue(y) = 1.0
```

To speed up computations during the lectures you may want to precompile the installed packages by running the following code:
```julia
using GLPKMathProgInterface
using Gadfly
using Interact
```

## More resources

We will not have the time to go through all of the basic syntax points of Julia. For more materials on learning Julia,
see [here](http://julialang.org/learning/). The JuMP documentation is located [here](http://www.juliaopt.org/JuMP.jl/0.14/).
During the lecture, we will be though through some of the examples in this repository and in these
[notebooks](http://nbviewer.jupyter.org/github/JuliaOpt/juliaopt-notebooks/tree/master/notebooks/).
