function jump_max_area_diameter(n)
    if mod(n,2)==0 && n>=6
        s = (2*sqrt(114)-7)/22;
        t = (3521*sqrt(114)-34010)/9196;
        c = (17328*(663157+3161*pi^2) - (1088031703 - 3918085*pi^2)*sqrt(114))/507398496;
        u = s*pi/n+t*pi/n^2-c*pi/n^3;
        v = (pi/2-u)/(n/2-1);
        d = asin(sin(u)+sin(u+3*v/2)/(2*cos(v/2)))-u-v;
        Al = sin(u)+sin(2*v)-sin(v+d)+(n/2-3)*(sin(v)-tan(v/2))+(cos(v-d)-cos(2*v)-1/2)*tan(v/2);
        Au = (n*sin(pi/n)-(n-1)*tan(pi/(2*n-2)))/2;
        a0 = zeros(Int(n/2));
        a0[1] = u;
        a0[2] = v+d;
        a0[3] = v-d;
        if n>=8
            a0[4:end] .= v;
        end
        model = Model(Ipopt.Optimizer)
        # Variables
        @variable(model, a[i = 1:Int(n/2)] >= 0, start = a0[i])
        # Area
        @NLobjective(model, Max, sin(a[1])+sum(sum((-1)^j*(sin(sum(a[k] for k in i-j:i+1))-
        sin(sum(a[k] for k in i-j:i))) for j in 0:i-2) for i = 2:Int(n/2-1)))
        # Optimal diameter graph
        @constraint(model, sum(a[k] for k in 1:Int(n/2)) == pi/2)
        @NLconstraint(model, sum((-1)^(j-1)*sin(sum(a[k] for k in 1:j)) for j in 1:Int(n/2-1)) == (-1)^(n/2)/2)
        @constraint(model, a[1] <= pi/6)
        @constraint(model, [i = 2:Int(n/2)], a[i] <= 2*a[1])
        # Bounds
        @NLconstraint(model, Al <= sin(a[1])+sum(sum((-1)^j*(sin(sum(a[k] for k in i-j:i+1))-
        sin(sum(a[k] for k in i-j:i))) for j in 0:i-2) for i = 2:Int(n/2-1)) <= Au)
        # Solving
        optimize!(model)
        A = objective_value(model)
        a = value.(a)
        ct = solve_time(model)
        return(A,a,ct)
    end
end