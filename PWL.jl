using JuMP, Cbc

# Generate supply, demand and piecewise linear objective
function generateData(snodes,dnodes,segments)

    supply     = rand(1:10,snodes)
    demand     = rand(1:10,dnodes)
    demandimbalance = sum(demand) - sum(supply)
    if demandimbalance > 0
        supply[snodes] += demandimbalance
    else
        demand[dnodes] -= demandimbalance
    end

    fvalues = [ Array{Float64}(segments+1) for i=1:snodes,j=1:dnodes ]
    dvalues = [ Array{Float64}(0) for i=1:snodes,j=1:dnodes ]
    for s in 1:snodes
        for d in 1:dnodes
            drange = min(supply[s],demand[d])
            delta = drange / (segments)
            dvalues[s,d] = [ (i-1)*delta for i=1:(segments+1) ]

            slopes = sort(rand(segments),rev=true)
            fvalues[s,d][1] = 0
            for i in 2:(segments+1)
                fvalues[s,d][i] = fvalues[s,d][i-1] + slopes[i-1]*(dvalues[s,d][i]-dvalues[s,d][i-1])
            end
        end
    end

    supply, demand, fvalues, dvalues
end

supply, demand, fvalues, dvalues = generateData(1,1,8)
snodes = length(supply); dnodes = length(demand);

function solveModel(supply, demand, fvalues, dvalues, pwlformulation,MIPsolver=CbcSolver())
    # Create JuMP model, using Cbc as the solver
    model = Model(solver=MIPsolver)

    # Define variables
    @variable(model, x[1:snodes,1:dnodes] >= 0)

    # Add constraints
    @constraints(model, begin
        supply_constraints[s=1:snodes], sum(x[s,d] for d in 1:dnodes) == supply[s]
        demand_constraints[d=1:dnodes], sum(x[s,d] for s in 1:snodes) == demand[d]
    end)

    # Define objective
    @objective(model, Min,
                           sum(piecewiselinear(model, x[s,d], dvalues[s,d], fvalues[s,d];method=pwlformulation)
                                                                             for s in 1:snodes, d in 1:dnodes))

    println("Solving for ",pwlformulation," ...")
    tic()
    status = solve(model)
    solvetime=toq();
    println("Solver status: ", status)
    println("Cost: ", getobjectivevalue(model))
    println("Solve time: ",solvetime)

end


using Gurobi


# Generate Data
supply, demand, fvalues, dvalues = generateData(10,10,32)
snodes = length(supply); dnodes = length(demand);

# Solve for four formulations using Gurobi
for formulation in [:CC,:MC,:SOS2,:Incremental,:Logarithmic,:DisaggLogarithmic,:ZigZag,:ZigZagInteger]
    solveModel(supply, demand, fvalues, dvalues, formulation,GurobiSolver(TimeLimit=10,OutputFlag=0))
end
