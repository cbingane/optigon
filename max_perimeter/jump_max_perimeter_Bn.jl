function jump_max_perimeter_Bn(n)
    if mod(log2(n),1)==0 && n>=8
        d = pi/n - asin(sin(2*pi/n)/2);
        Ll = 2*n*sin(pi/(2*n))*cos(d/2);
        Lu = 2*n*sin(pi/(2*n));
        #model = Model(Ipopt.Optimizer)
        #model = Model(() -> AmplNLWriter.Optimizer(Couenne_jll.amplexe))
        model = Model(() -> AmplNLWriter.Optimizer("C:/Users/chris/couenne-win64/couenne.exe"))
        # Variables
        @variable(model, a[i = 1:Int(n/4+1)] >= 0, start = pi/n+(-1)^(i-1)*d)
        # Perimeter
        @NLobjective(model, Max, 4*sin(a[1]/2)+sum(8*sin(a[i]/2) for i in 2:Int(n/4))+4*sin(a[Int(n/4+1)]/2))
        # Diameter graph
        @constraint(model, a[1] + sum(2*a[k] for k in 2:Int(n/4)) + a[Int(n/4+1)] == pi/2)
        @NLconstraint(model, sin(a[1])+sum(cos((i-1)*pi)*sin(a[1]+
        sum(2*a[j] for j in 2:i)) for i in 2:Int(n/4)) == -1/2)
        @constraint(model, [i = 1:Int(n/4)], a[i] <= pi/6)
        @constraint(model, a[Int(n/4+1)] <= pi/3)
        # Bounds
        @NLconstraint(model, Ll <= 4*sin(a[1]/2)+sum(8*sin(a[i]/2) for i in 2:Int(n/4))+4*sin(a[Int(n/4+1)]/2) <= Lu)
        # Solving
        optimize!(model)
        L = objective_value(model)
        a = value.(a)
        ct = solve_time(model)
        return(L,a,ct)
    end
end