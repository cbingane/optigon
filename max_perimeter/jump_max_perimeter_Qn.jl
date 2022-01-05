function jump_max_area_perimeter_Qn(n)
    if mod(log2(n),1)==0 && n>=4
        d = pi/4-asin(cos(pi/n)/sqrt(2));
        Ll = 2*n*sin(pi/(2*n))*cos(d/2);
        Lu = 2*n*sin(pi/(2*n));
        model = Model(Ipopt.Optimizer)
        #model = Model(() -> AmplNLWriter.Optimizer(Couenne_jll.amplexe))
        # Variables
        @variable(model, a[i = 1:Int(n/2)] >= 0, start = pi/n+(-1)^i*d)
        # Perimeter
        @NLobjective(model, Max, sum(4*sin(a[i]/2) for i in 1:Int(n/2)))
        # Optimal diameter graph
        @constraint(model, sum(a[i] for i in 1:Int(n/2)) == pi/2)
        @NLconstraint(model, sum((-1)^(i-1)*sin(sum(a[j] for j in 1:i)) for i in 1:Int(n/2-1)) == 1/2)
        @constraint(model, a[1] <= pi/6)
        @constraint(model, [i = 2:Int(n/2)], a[i] <= pi/3)
        # Bounds
        @NLconstraint(model, Ll <= sum(4*sin(a[i]/2) for i in 1:Int(n/2)) <= Lu)
        # Solving
        optimize!(model)
        L = objective_value(model)
        a = value.(a)
        ct = solve_time(model)
        return(L,a,ct)
    end
end