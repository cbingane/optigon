# 'cstrt_bwinja_ngon' provides the vertices coordinates (a,b) of a small
# n-gon for n = 2m and m >= 3
function cstrt_bwinja_ngon(n)
if mod(n,2) == 0 && n>=6
    # initialization
    v(u) = (pi/2-u)/(n/2-1);
    d(u) = asin(sin(u)+sin(u+3*v(u)/2)/(2*cos(v(u)/2)))-u-v(u);
    A(u) = sin(u)+sin(2*v(u))-sin(v(u)+d(u))+
    (n/2-3)*(sin(v(u))-tan(v(u)/2))+(cos(v(u)-d(u))-cos(2*v(u))-1/2)*tan(v(u)/2);
    F(u) = -A(u);
    res = optimize(F,pi/(2*n-2),pi/n);
    u0 = Optim.minimizer(res);
    v0 = (pi/2-u0)/(n/2-1);
    d0 = asin(sin(u0)+sin(u0+3*v0/2)/(2*cos(v0/2)))-u0-v0;
    x = zeros(Int(n/2));
    x[1] = u0;
    x[2] = v0+d0;
    x[3] = v0-d0;
    if n >= 8
        x[4:end] .= v0;
    end
    # construction
    a = zeros(n-1); b = zeros(n-1);
    a[Int(n/2-1)] = sin(x[1]); a[Int(n/2+1)] = -a[Int(n/2-1)];
    b[Int(n/2-1)] = cos(x[1]); b[Int(n/2+1)] = b[Int(n/2-1)];
    for i = 2:Int(n/2-1)
        a[Int(mod(i*(n-1)/2-mod(i,2)/2,n))] = a[Int(mod((i-1)*(n-1)/2-mod(i-1,2)/2,n))]-(-1)^i*sin(sum(x[1:i]));
        a[Int(n-mod(i*(n-1)/2-mod(i,2)/2,n))] = -a[Int(mod(i*(n-1)/2-mod(i,2)/2,n))];
        b[Int(mod(i*(n-1)/2-mod(i,2)/2,n))] = b[Int(mod((i-1)*(n-1)/2-mod(i-1,2)/2,n))]-(-1)^i*cos(sum(x[1:i]));
        b[Int(n-mod(i*(n-1)/2-mod(i,2)/2,n))] = b[Int(mod(i*(n-1)/2-mod(i,2)/2,n))];
    end
    a[Int(n/2)] = 0; b[Int(n/2)] = 1;
    return(a,b)
end
end