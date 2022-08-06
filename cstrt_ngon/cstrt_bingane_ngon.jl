# 'cstrt_bingane_ngon' provides the vertices coordinates (a,b) of a convex
# small n-gon Bn for n = 2^s and s >= 3
# The 8-gon has maximal width among all convex small 8-gons
# Please, cite...
#   C. Bingane. Tight bounds on the maximal perimeter and the maximal width
#   of convex small polygons. Journal of Global Optimization, 2022.
function cstrt_bingane_ngon(n)
    if mod(log2(n),1)==0 && n>=8
        # initialization
        Pi = big(pi);
        d = Pi/n - asin(sin(2*Pi/n)/2);
        # construction
        a = Vector{BigFloat}(undef,n-1); b = Vector{BigFloat}(undef,n-1);
        for i = 1:Int(n/4)
            a[Int(mod(i*(n/2-1),n))] = sin(2*i*Pi/n)*sin(d + (-1)^(i-1)*pi/n)/sin(2*Pi/n);
            a[Int(n-mod(i*(n/2-1),n))] = -a[Int(mod(i*(n/2-1),n))];
            b[Int(mod(i*(n/2-1),n))] = (sin(Pi/n - d) + cos(2*i*Pi/n)*sin(d + (-1)^(i-1)*Pi/n))/sin(2*Pi/n);
            b[Int(n-mod(i*(n/2-1),n))] = b[Int(mod(i*(n/2-1),n))];
        end
        for i = 1:Int(n/4-1)
            a[Int(mod(i*(n/2-1)+n/2,n))] = a[Int(mod(i*(n/2-1),n))] + (-1)^i*sin(2*i*Pi/n);
            a[Int(n-mod(i*(n/2-1)+n/2,n))] = -a[Int(mod(i*(n/2-1)+n/2,n))];
            b[Int(mod(i*(n/2-1)+n/2,n))] = b[Int(mod(i*(n/2-1),n))] + (-1)^i*cos(2*i*Pi/n);
            b[Int(n-mod(i*(n/2-1)+n/2,n))] = b[Int(mod(i*(n/2-1)+n/2,n))];
        end
        a[Int(n/2)] = 0; b[Int(n/2)] = 1;
        return(a,b)
    else
        throw(DomainError(n,"n must be a power of 2 >= 8"))
    end
end