# 'jump_max_area_xy_rstr_quad' solves a convex restriction of the maximal
# area problem constructed around an n-gon (a,b)
function jump_max_area_xy_rstr_quad(n,a,b)
    if mod(n,2)==0 && n>=6
        #model = Model(Ipopt.Optimizer)
        #model = Model(Mosek.Optimizer)
        model = Model(with_optimizer(Mosek.Optimizer,MSK_IPAR_INTPNT_SOLVE_FORM=MSK_SOLVE_DUAL))
        # Variables
        @variable(model, x[i = 1:n-1], start = a[i])
        @variable(model, y[i = 1:n-1] >= 0, start = b[i])
        @variable(model, u[i = 1:n-2] >= 0, start = (a[i]*b[i+1]-b[i]*a[i+1])/2)
        # Area
        @objective(model, Max, sum(u[i]  for i in 1:n-2))
        # Unit-diameter
        @constraint(model, [i = 1:n-1], x[i]^2+y[i]^2 <= 1)
        @constraint(model, [i = 2:n-1, j = 1:i-1], (x[i]-x[j])^2+(y[i]-y[j])^2 <= 1)
        # Counterclockwise order
        @constraint(model, [i = 1:n-2], (x[i]-y[i+1])^2+(y[i]+x[i+1])^2+8*u[i] <=
        2*(a[i]+b[i+1])*(x[i]+y[i+1])-(a[i]+b[i+1])^2+2*(b[i]-a[i+1])*(y[i]-x[i+1])-(b[i]-a[i+1])^2)
        # Solving
        optimize!(model)
        x = value.(x)
        y = value.(y)
        ct = solve_time(model)
        return(x,y,ct)
    end
end