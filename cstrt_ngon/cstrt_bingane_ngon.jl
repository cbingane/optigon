# 'cstrt_bingane_ngon' provides the vertices coordinates (a,b) of a convex
# small n-gon Bn for n = 2^s and s >= 3
# The 8-gon has maximal width among all convex small 8-gons
# Please, cite...
#   C. Bingane. Tight bounds on the maximal perimeter and the maximal width
#   of convex small polygons. Journal of Global Optimization, 2022.
function cstrt_bingane_ngon(n)
if mod(log2(n),1)==0 && n>=8
    # initialization
    d = pi/n - asin(sin(2*pi/n)/2);
    # construction
    a = zeros(n-1); b = zeros(n-1);
    for i = 1:Int(n/4)
        a[Int(mod(i*(n/2-1),n))] = sin(2*i*pi/n)*sin(d + (-1)^(i-1)*pi/n)/sin(2*pi/n);
        a[Int(n-mod(i*(n/2-1),n))] = -a[Int(mod(i*(n/2-1),n))];
        b[Int(mod(i*(n/2-1),n))] = (sin(pi/n - d) + cos(2*i*pi/n)*sin(d + (-1)^(i-1)*pi/n))/sin(2*pi/n);
        b[Int(n-mod(i*(n/2-1),n))] = b[Int(mod(i*(n/2-1),n))];
    end
    for i = 1:Int(n/4-1)
        a[Int(mod(i*(n/2-1)+n/2,n))] = a[Int(mod(i*(n/2-1),n))] + (-1)^i*sin(2*i*pi/n);
        a[Int(n-mod(i*(n/2-1)+n/2,n))] = -a[Int(mod(i*(n/2-1)+n/2,n))];
        b[Int(mod(i*(n/2-1)+n/2,n))] = b[Int(mod(i*(n/2-1),n))] + (-1)^i*cos(2*i*pi/n);
        b[Int(n-mod(i*(n/2-1)+n/2,n))] = b[Int(mod(i*(n/2-1)+n/2,n))];
    end
    a[Int(n/2)] = 0; b[Int(n/2)] = 1;
    return(a,b)
end
end