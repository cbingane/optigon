# 'jump_max_width_F8' computes the maximal width of an equilateral
# small 8-gon
# Please, cite...
#   C. Bingane and C. Audet. The equilateral small octagon of maximal
#   width. Mathematics of Computation, 2022.

#model = Model(Ipopt.Optimizer)
#model = Model(() -> AmplNLWriter.Optimizer(Couenne_jll.amplexe))
model = Model(() -> AmplNLWriter.Optimizer("C:/Users/chris/couenne-win64/couenne.exe"))
n = 8;
# Optimal height graph
k = [5,5,6,0,0,2,3,3];
# Variables
@variable(model, x[i = 1:n-1])
@variable(model, y[i = 1:n-1] >= 0)
@variable(model, 0.384462 <= c <= 0.386952)
@variable(model, h[i = 1:n])
@variable(model, W >= 0.953776)
# Width
@objective(model, Max, W)
@constraint(model, [i = 1:n], h[i] >= W)
# Diameter
@NLconstraint(model, [i = 1:n-1], x[i]^2+y[i]^2 <= 1)
@NLconstraint(model, [i = 2:n-1, j = 1:i-1], (x[i]-x[j])^2+(y[i]-y[j])^2 <= 1)
# Counterclockwise order
@NLconstraint(model, [i = 1:n-2], x[i]*y[i+1] >= y[i]*x[i+1])
# Sides
@NLconstraint(model, x[1]^2+y[1]^2 == c^2)
@NLconstraint(model, [i = 2:n-1], (x[i]-x[i-1])^2+(y[i]-y[i-1])^2 == c^2)
@NLconstraint(model, x[n-1]^2+y[n-1]^2 == c^2)
# Heights
@NLconstraint(model, c*h[1] == -x[k[1]]*(y[1]-y[k[1]])+y[k[1]]*(x[1]-x[k[1]]))
@NLconstraint(model, [i = 2:n-1; k[i] == 0], c*h[i] == x[i-1]*y[i]-y[i-1]*x[i])
@NLconstraint(model, [i = 2:n-1; k[i] != 0], c*h[i] == (x[i-1]-x[k[i]])*(y[i]-y[k[i]])-(y[i-1]-y[k[i]])*(x[i]-x[k[i]]))
@NLconstraint(model, c*h[n] == -(x[n-1]-x[k[n]])*y[k[n]]+(y[n-1]-y[k[n]])*x[k[n]])
# Break symmetry
@constraint(model, x[Int(n/2)] == 0)
# Solving
optimize!(model)
W = objective_value(model);
h = value.(h);
c = value.(c);
x = value.(x);
y = value.(y);
ct = solve_time(model);