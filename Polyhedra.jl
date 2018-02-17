using JuMP, Polyhedra, CDDLib

m = Model()
n = 3
@variable(m, x[1:n])
@variable(m, y[1:n] >= 0)
@constraint(m,[i=1:n],  x[i] <= y[i])
@constraint(m,[i=1:n], -x[i] <= y[i])
@constraint(m, sum(y) <= 1)
m

poly = polyhedron(m, CDDLibrary(:exact))

poly_x = eliminate(poly, n+1:2n)

@show SimpleHRepresentation(poly_x)
removehredundancy!(poly_x)
@show SimpleHRepresentation(poly_x);

@show SimpleVRepresentation(poly_x)
removevredundancy!(poly_x)
@show SimpleVRepresentation(poly_x);
