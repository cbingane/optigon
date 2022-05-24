# 'cstrt_hansen_ngon' provides the vertices coordinates (a,b) of a small
# n-gon for n = 2^s and s >= 2
# The 4-gon has maximal width among all 4-gons with fixed perimeter
function cstrt_hansen_ngon(n)
if mod(log2(n),1)==0 && n>=4
    # initialization
    v(u) = (pi/2-u)/(n/2-1);
    W(u) = 2*cos(2*u-v(u))*cos(u)*cos(v(u)-u)/(n*sin(2*u)+(n-2)*sin(2*v(u)-2*u));
    F(u) = -W(u);
    res = Optim.optimize(F,pi/(2*n-2),pi/n);
    u0 = Optim.minimizer(res);
    v0 = (pi/2-u0)/(n/2-1);
    d0 = 2*cos(2*u0-v0)*cos(u0)/(n*sin(2*u0)+(n-2)*sin(2*v0-2*u0));
    D0 = d0*cos(v0-u0)/cos(u0);
    z0 = d0*cos(v0-u0)/cos(2*u0-v0);
    r = ones(Int(n/2+1));
    for i = 1:Int(n/4)
        r[2*i-1] = d0/D0;
    end
    r[Int(n/2+1)] = z0/D0;
    x = zeros(Int(n/2));
    x[1] = u0;
    x[2:end] .= v0;
    # construction
    a = zeros(n-1); b = zeros(n-1);
    a[Int(n/2-1)] = r[1]*sin(x[1]); a[Int(n/2+1)] = -a[Int(n/2-1)];
    b[Int(n/2-1)] = r[1]*cos(x[1]); b[Int(n/2+1)] = b[Int(n/2-1)];
    for i = 2:Int(n/2-1)
        a[Int(mod(i*(n-1)/2-mod(i,2)/2,n))] = a[Int(mod((i-1)*(n-1)/2-mod(i-1,2)/2,n))] - (-1)^i*r[i]*sin(sum(x[1:i]));
        a[Int(n-mod(i*(n-1)/2-mod(i,2)/2,n))] = -a[Int(mod(i*(n-1)/2-mod(i,2)/2,n))];
        b[Int(mod(i*(n-1)/2-mod(i,2)/2,n))] = b[Int(mod((i-1)*(n-1)/2-mod(i-1,2)/2,n))] - (-1)^i*r[i]*cos(sum(x[1:i]));
        b[Int(n-mod(i*(n-1)/2-mod(i,2)/2,n))] = b[Int(mod(i*(n-1)/2-mod(i,2)/2,n))];
    end
    a[Int(n/2)] = 0; b[Int(n/2)] = r[Int(n/2+1)];
    return(a,b)
end
end