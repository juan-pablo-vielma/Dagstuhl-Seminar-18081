using JuMP, PolyJuMP, MultivariatePolynomials, DynamicPolynomials, SumOfSquares

type Box
    xl::Float64
    xu::Float64
    yl::Float64
    yu::Float64
end

function solveModel(npieces,degree,boxes::Array{Box},domain::Box,startcondition,endcondition,misdpsolver,tol_int = 1e-9, tol_feas = 1e-9, tol_gap = 1e-6)

    model = SOSModel(solver=misdpsolver)

    tol_nonneg = tol_feas # Tolerance for polynomial nonnegativity

    (X₀, X₀′, X₀′′) = startcondition
    (X₁, X₁′, X₁′′) = endcondition

    # Discretize time into npieces+1 times
    Tmin = 0.
    Tmax = 1.
    T = linspace(Tmin, Tmax, npieces+1)

    # Polynomials are a function of t
    @polyvar(t)
    Z = monomials(t, 0:degree)

    # Binary variables to choose safe regions
    # variables are indexed by trajectory piece (integer 1:npieces) and boxes (of type Box)
    @variable(model, H[1:npieces,boxes], Bin)

    # Big-M values
    (Mxl, Mxu, Myl, Myu) = (domain.xl, domain.xu, domain.yl, domain.yu)
    p = Dict()
    for j in 1:npieces
        @constraint(model, sum(H[j,box] for box in boxes) == 1)

        # Polynomial variables = a variable for the coefficients of each monomial in Z
        p[(:x,j)] = @variable(model, _, PolyJuMP.Poly(Z), basename="px$j")
        p[(:y,j)] = @variable(model, _, PolyJuMP.Poly(Z), basename="py$j")

        # Constraints to choose safe region
        for box in boxes
            xl, xu, yl, yu = box.xl, box.xu, box.yl, box.yu
            @assert xl >= Mxl
            @constraint(model, p[(:x,j)] >= Mxl + (xl-Mxl)*H[j,box], domain = (@set t >= T[j] && t <= T[j+1]))
            @assert xu <= Mxu
            @constraint(model, p[(:x,j)] <= Mxu + (xu-Mxu)*H[j,box], domain = (@set t >= T[j] && t <= T[j+1]))
            @assert yl >= Myl
            @constraint(model, p[(:y,j)] >= Myl + (yl-Myl)*H[j,box], domain = (@set t >= T[j] && t <= T[j+1]))
            @assert yu <= Myu
            @constraint(model, p[(:y,j)] <= Myu + (yu-Myu)*H[j,box], domain = (@set t >= T[j] && t <= T[j+1]))
        end
    end

    # Boundary and interstitial smoothing conditions
    for axis in (:x,:y)
        @constraint(model,               p[(axis,1)       ](t=>Tmin) == X₀[axis])
        @constraint(model, differentiate(p[(axis,1)], t   )(t=>Tmin) == X₀′[axis])
        @constraint(model, differentiate(p[(axis,1)], t, 2)(t=>Tmin) == X₀′′[axis])

        for j in 1:npieces-1
            @constraint(model,               p[(axis,j)       ](t=>T[j+1]) ==               p[(axis,j+1)       ](t=>T[j+1]))
            @constraint(model, differentiate(p[(axis,j)], t   )(t=>T[j+1]) == differentiate(p[(axis,j+1)], t   )(t=>T[j+1]))
            @constraint(model, differentiate(p[(axis,j)], t, 2)(t=>T[j+1]) == differentiate(p[(axis,j+1)], t, 2)(t=>T[j+1]))
        end

        @constraint(model, p[(axis,npieces)](t=>Tmax) == X₁[axis])
    end

    # Objective function
    @variable(model, γ[keys(p)] ≥ 0)
    for (key,val) in p
        @constraint(model, γ[key] ≥ norm(differentiate(val, t, 3)))
    end
    @objective(model, Min, sum(γ))

    solve(model)

    # create and return a function that evaluates the piecewise polynomial anwer
    PP = Dict(key => getvalue(p[key]) for key in keys(p))
    HH = getvalue(H)
    function eval_poly(r)
        for i in 1:npieces
            if T[i] <= r <= T[i+1]
                return PP[(:x,i)](t=>r), PP[(:y,i)](t=>r)
                break
            end
        end
        error("Time $r out of interval [$(minimum(T)),$(maximum(T))]")
    end
end

# Define instance
npieces = 8 # Number of trajectory pieces
degree = 3 # degree of polynomial trajectories
boxes = Box[
    Box(0.0,1.0,0.0,0.3),
    Box(0.8,1.7,0.1,0.3),
    Box(1.4,1.9,0.2,0.4),
    Box(1.0,1.7,0.3,0.5),
    Box(0.5,1.4,0.4,0.6),
    Box(0.0,1.0,0.5,0.7),
    Box(0.2,1.0,0.6,0.8),
    Box(0.5,1.3,0.7,0.9),
    Box(1.0,2.0,0.7,1.0),
    ]
domain = Box(0,2,0,1)
# Start location, velocity and acceleration
X₀   = Dict(:x=>0, :y=>0)
X₀′  = Dict(:x=>1, :y=>0)
X₀′′ = Dict(:x=>0, :y=>0)
# End location, velocity and acceleration
X₁   = Dict(:x=>2, :y=>1)
X₁′  = Dict(:x=>1, :y=>0)
X₁′′ = Dict(:x=>0, :y=>0)

# Setup Pajarito MISDP solver
using Pajarito, CPLEX, Mosek
MIPsolver = CplexSolver()
SDPsolver = MosekSolver(MSK_IPAR_LOG=0,MSK_DPAR_INTPNT_CO_TOL_PFEAS=1e-6,MSK_DPAR_INTPNT_CO_TOL_DFEAS=1e-6,
    MSK_DPAR_INTPNT_CO_TOL_REL_GAP=1e-5,MSK_DPAR_INTPNT_TOL_INFEAS=1e-8,MSK_DPAR_INTPNT_CO_TOL_MU_RED=1e-6)
misdpsolver = PajaritoSolver(mip_solver=MIPsolver, mip_subopt_solver=MIPsolver, cont_solver=SDPsolver, mip_solver_drives=true, cut_zero_tol=1e-6, log_level=3, rel_gap=1e-3, solve_relax=false)

# Solve
eval_poly = solveModel(npieces,degree,boxes,domain,(X₀, X₀′, X₀′′),(X₁, X₁′, X₁′′),misdpsolver)

eval_poly(0.4)
