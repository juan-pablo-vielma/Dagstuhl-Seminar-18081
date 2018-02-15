#  JuMP Demo at Dagstuhl Seminar 18081

This site contains materials for a demonstration of [JuMP](https://github.com/JuliaOpt/JuMP.jl) at the [Dagstuhl Seminar 18081: Designing and Implementing Algorithms for Mixed-Integer Nonlinear Optimization](http://www.dagstuhl.de/en/program/calendar/semhp/?semnr=18081)

## Install Julia

You should use the latest version of Julia v0.6.2.
Binaries of Julia for all platforms are available [here](http://julialang.org/downloads/).


## Install Basic JuMP packages

To install [JuMP](https://github.com/JuliaOpt/JuMP.jl) and the open-source LP/MIP solvers [Clp](https://projects.coin-or.org/Clp) and [Cbc](https://projects.coin-or.org/Cbc)  run the following code:
```julia
Pkg.add("JuMP")
Pkg.add("Clp")
```
If you have a previous installation of Julia,
be sure to update your packages to the latest version by running ``Pkg.update()``.

To test that your installation is working, run the following code (the first time you run the code you may see the message like "INFO: Precompiling stale cache ..." for a few seconds):

```julia
using JuMP, Cbc
m = Model(solver=CbcSolver())
@variable(m, x >= 0, Int)
@variable(m, y >= 0)
@constraint(m, 2x + y <= 1)
@objective(m, Max, x+y)
status = solve(m)
@show status; @show getvalue(y);
```

The output should be:

```
status = :Optimal
getvalue(y) = 1.0
```

## Install Visualization Package

For visualization we will use the [Plots](https://github.com/JuliaPlots/Plots.jl) package, which can be installed by running:
```julia
Pkg.add("Plots")
```

## Install Visualization Package

One of our demonstrations will use advanced mixed integer programming (MIP) formulations for piecewise linear functions. These can be easily accessed through the [PiecewiseLinearOpt](https://github.com/joehuchette/PiecewiseLinearOpt.jl) package, which can be installed by running:
```julia
Pkg.add("PiecewiseLinearOpt")
```

## Install Polyhedra Package

Another demonstration will use the package [Polyhedra](https://github.com/JuliaPolyhedra/Polyhedra.jl), which can be used to obtain the extreme points of the feasible region of a linear programming JuMP model by calling a library such as [CDD](https://www.inf.ethz.ch/personal/fukudak/cdd_home/). The Polyhedra package and the interphase for CDD can be installed by running:

```julia
Pkg.add("Polyhedra")
Pkg.add("CDDLib")
```

Package CDDLib downloads binaries for Windows and compiles CDDLib from sources for Linux and Mac. Install from sources just requires the presence of a C compiler and Linux additionally requires the GMP header files (see details in the [CDDLib page](https://github.com/JuliaPolyhedra/CDDLib.jl)).

## Install Packages for Handling Polynomials

For our advanced demo we will solve mixed integer optimal control problems where the variables are polynomials. To solve these problems we construct a finite dimensional approximation using sum-of-square (SOS) programming that can be formulated as a mixed integer semidefinite programming (MISDP) problem. This requires several packages to handle polynomials ([MultivariatePolynomials](https://github.com/JuliaAlgebra/MultivariatePolynomials.jl) and [DynamicPolynomials](https://github.com/JuliaAlgebra/DynamicPolynomials.jl)), extend JuMP to consider polynomial variables ([PolyJuMP](https://github.com/JuliaOpt/PolyJuMP.jl)), and construct the SOS approximations ([SumOfSquares](https://github.com/JuliaOpt/SumOfSquares.jl)).

To install all these packages run the following code:
```julia
Pkg.add("PolyJuMP")
Pkg.add("MultivariatePolynomials")
Pkg.add("DynamicPolynomials")
Pkg.add("SumOfSquares")
```

To speed up computations during the lectures you may want to precompile the installed packages by running the following code:
```julia
using PolyJuMP, MultivariatePolynomials, DynamicPolynomials, SumOfSquares
```

## Install MICP Solver Pajarito

To solve the MISDP formulation of the control problem we will use the Julia-written mixed integer convex programming (MICP) solver [Pajarito](https://github.com/JuliaOpt/Pajarito.jl).

To install Pajarito run the following code:
```julia
Pkg.add("Pajarito")
```

Pajarito requires the presence  of a linear mixed integer programming (MIP) solver and a semidefinite programming (SDP). While Pajarito can use Cbc as a MIP solver, we recommend also installing a commercial MIP solver.  We include instructions to install an SDP solver and a commercial MIP solver in the following sections.


## Install SDP Solver

We recommend [Mosek](https://www.mosek.com) v8 as a SDP solver.

To install Mosek in Windows or Linux run the following code:
```julia
Pkg.add("Mosek")
```
To install Mosek in Mac run the following code (there is a small [bug-fix](https://github.com/JuliaOpt/Mosek.jl/issues/134) not yet in the release version):
```julia
Pkg.add("Mosek")
Pkg.checkout("Mosek","b0.8")
Pkg.build("Mosek")
```
In both cases the package install automatically downloads the Mosek binaries, but you need to request an academic license file using this [form](https://license.mosek.com/academic/).

If you do not have access to a Mosek license we recommend installing the open-source SDP solver [SCS](https://github.com/cvxgrp/scs) by running the following code:
```julia
Pkg.add("SCS")
```
This package install also automatically downloads the SCS source and compiles the associated binaries.

You may also want to install [ECOS](https://github.com/embotech/ecos) to solve SOCP problems. For this run the following code:
```julia
Pkg.add("ECOS")
```
This package install also automatically downloads the ECOS source and compiles the associated binaries.


## Install Commercial MIP Solver

We also recommend installing a commercial MIP solver. The install instructions below assume you already have the corresponding binaries and license files installed in your system.

To make [CPLEX](https://www.ibm.com/analytics/data-science/prescriptive-analytics/cplex-optimizer) available to JuMP run the following code:
```julia
Pkg.add("CPLEX")
```

To make [Gurobi](http://www.gurobi.com) available to JuMP run the following code:
```julia
Pkg.add("Gurobi")
```

## More resources

We will not have the time to go through all of the basic syntax points of Julia. For more materials on learning Julia,
see [here](http://julialang.org/learning/). The JuMP documentation is located [here](http://www.juliaopt.org/JuMP.jl/0.18/).
