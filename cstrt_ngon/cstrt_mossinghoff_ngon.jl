# 'cstrt_mossinghoff_ngon' provides the vertices coordinates (a,b) of a small
# n-gon for n = 2m and m >= 3
# Please, cite...
#   C. Bingane. Tight bounds on the maximal area of small polygons:
#   Improved Mossinghoff polygons. Discrete & Computational Geometry, 2022.
#   M. J. Mossinghoff. Isodiametric problems for polygons. Discrete &
#   Computational Geometry, 36(2): 363-379, 2006.
function cstrt_mossinghoff_ngon(n)
if mod(n,2) == 0 && n>=6
    # initialization
    s = (2*sqrt(114)-7)/22;
    t = (3521*sqrt(114)-34010)/9196;
    c = (17328*(663157+3161*pi^2) - (1088031703 - 3918085*pi^2)*sqrt(114))/507398496;
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