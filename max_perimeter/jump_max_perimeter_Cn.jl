function jump_max_perimeter_Dn(n)
    if mod(log2(n),1)==0 && n>=16
        d = atan(tan(2*pi/n)*tan(pi/n))-asin(sin(2*pi/n)^2/sqrt(2+2*cos(2*pi/n)*cos(4*pi/n)));
        Ll = 2*n*sin(pi/(2*n))*cos(d/2);
        Lu = 2*n*sin(pi/(2*n));
        c = zeros(Int(3*n/8));
        for i = 1:Int(3*n/8)
            c[i] = 1 + round(mod(i,3)/3);
        end
        #model = Model(Ipopt.Optimizer)
        #model = Model(() -> AmplNLWriter.Optimizer(Couenne_jll.amplexe))
        model = Model(() -> AmplNLWriter.Optimizer("C:/Users/chris/couenne-win64/couenne.exe"))
        # Variables
        @variable(model, a[i = 1:Int(3*n/8)] >= 0, start = pi/n+(-1)^(i-1)*d)
        # Perimeter
        @NLobjective(model, Max, sum(4*c[i]*sin(a[i]/2) for i in 1:Int(3*n/8)))
        # Diameter graph
        @constraint(model, sum(c[i]*a[i] for i in 1:Int(3*n/8)) == pi/2)
        @NLconstraint(model, sum(cos((i-1)*pi)*sin(sum(c[j]*a[j] for j in 1:i)) for i in 1:Int(3*n/8-1)) == 1/2)
        @constraint(model, [i = 1:Int(3*n/8)], c[i]*a[i] <= pi/3)
        # Bounds
        @NLconstraint(model, Ll <= sum(4*c[i]*sin(a[i]/2) for i in 1:Int(3*n/8)) <= Lu)
        # Solving
        optimize!(model)
        L = objective_value(model)
        a = value.(a)
        ct = solve_time(model)
        return(L,a,ct)
    end
end
