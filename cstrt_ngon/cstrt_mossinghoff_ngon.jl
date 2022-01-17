# 'cstrt_bieri_ngon' provides the vertices coordinates (a,b) of a small
# n-gon for n = 2m and m >= 3
function cstrt_bieri_ngon(n)
if mod(n,2) == 0 && n>=6
    # initialization
    s = (2*sqrt(114)-7)/22;
    t = (84*s^2-272*s+175)/(4*(22*s+7));
    c = (7792*s^4+16096*s^3+2568*s^2-6248*s+223)*pi^2/(768*(22*s+7))-(226*s^2+84*s*t-22*t^2-542*s-136*t+303)/(2*(22*s+7));
    u = s*pi/n + t*pi/n^2 - c*pi/n^3;
    v = (pi/2-u)/(n/2-1);
    d = asin(sin(u)+sin(u+3*v/2)/(2*cos(v/2)))-u-v;
    x = zeros(Int(n/2));
    x[1] = u;
    x[2] = v+d;
    x[3] = v-d;
    if n >= 8
        x[4:end] .= v;
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