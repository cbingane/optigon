# 'cstrt_messine_ngon' provides the vertices coordinates (a,b) of a convex
# small n-gon for n = 2^s and s >= 2
# The 4-gon has maximal area, maximal perimeter and maximal width among all
# convex small 4-gons
# Please, cite...
#   C. Bingane. Tight bounds on the maximal perimeter and the maximal width
#   of convex small polygons. Journal of Global Optimization, 2022.
function cstrt_messine_ngon(n)
if mod(log2(n),1)==0 && n>=4
    # initialization
    d = pi/4 - asin(cos(pi/n)/sqrt(2));
    x = zeros(Int(n/2));
    for i = 1:Int(n/2)
        x[i] = pi/n + (-1)^i*d;
    end
    # construction
    a = zeros(n-1); b = zeros(n-1);
    a[Int(n/2-1)] = sin(x[1]); a[Int(n/2+1)] = -a[Int(n/2-1)];
    b[Int(n/2-1)] = cos(x[1]); b[Int(n/2+1)] = b[Int(n/2-1)];
    for i = 2:Int(n/2-1)
        a[Int(mod(i*(n-1)/2-mod(i,2)/2,n))] = a[Int(mod((i-1)*(n-1)/2-mod(i-1,2)/2,n))] - (-1)^i*sin(sum(x[1:i]));
        a[Int(n-mod(i*(n-1)/2-mod(i,2)/2,n))] = -a[Int(mod(i*(n-1)/2-mod(i,2)/2,n))];
        b[Int(mod(i*(n-1)/2-mod(i,2)/2,n))] = b[Int(mod((i-1)*(n-1)/2-mod(i-1,2)/2,n))] - (-1)^i*cos(sum(x[1:i]));
        b[Int(n-mod(i*(n-1)/2-mod(i,2)/2,n))] = b[Int(mod(i*(n-1)/2-mod(i,2)/2,n))];
    end
    a[Int(n/2)] = 0; b[Int(n/2)] = 1;
    return(a,b)
end
end