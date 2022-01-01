# 'cstrt_xiong_ngon' provides the vertices coordinates (a,b) of a small
# n-gon for n = 2m and m >= 4
function cstrt_xiong_ngon(n)
if mod(n,2) == 0 && n>=8
    # initialization
    v(u) = (pi/2-u[1]-2*u[2])/(n/2-3);
    d(u) = asin(sin(u[1])+sin(u[1]+2*u[2]-v(u)/2)/(2*cos(v(u)/2)))-u[1]-u[2];
    A(u) = sin(u[1])+sin(2*u[2])-sin(u[2]+d(u))+
    (n/2-3)*(sin(v(u))-tan(v(u)/2))+(cos(u[2]-d(u))-cos(2*u[2])-1/2)*tan(v(u)/2);
    F(u) = -A(u);
    res = optimize(F,[pi/(2*n-2), pi/(n-1)],LBFGS());
    u0, w0 = Optim.minimizer(res);
    v0 = (pi/2-u0-2*w0)/(n/2-3);
    d0 = asin(sin(u0)+sin(u0+2*w0-v0/2)/(2*cos(v0/2)))-u0-w0;
    x = zeros(Int(n/2));
    x[1] = u0;
    x[2] = w0+d0;
    x[3] = w0-d0;
    x[4:end] .= v0;
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