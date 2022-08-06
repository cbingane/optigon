# 'cstrt_taylor_ngon' provides the vertices coordinates (a,b) of a convex
# small n-gon for n = 2^s and s >= 2
# The 4-gon has maximal area, maximal perimeter and maximal width among all
# convex small 4-gons
# Please, cite...
#   C. Bingane. Tight bounds on the maximal perimeter and the maximal width
#   of convex small polygons. Journal of Global Optimization, 2022.
function cstrt_taylor_ngon(n)
    if mod(log2(n),1)==0 && n>=4
        # initialization
        m = div(n,2);
        Pi = big(pi);
        d = Pi/4 - asin(cos(Pi/n)/sqrt(big(2)));
        x = Vector{BigFloat}(undef,m);
        for i = 1:m
            x[i] = Pi/n + (-1)^i*d;
        end
        # construction
        a = Vector{BigFloat}(undef,n-1); b = Vector{BigFloat}(undef,n-1);
        a[m-1] = sin(x[1]); a[m+1] = -a[m-1];
        b[m-1] = cos(x[1]); b[m+1] = b[m-1];
        if n>4
            for i = 2:m-1
                a[mod(div(i*(n-1)-mod(i,2),2),n)] = a[mod(div((i-1)*(n-1)-mod(i-1,2),2),n)] - (-1)^i*sin(sum(x[1:i]));
                a[n-mod(div(i*(n-1)-mod(i,2),2),n)] = -a[mod(div(i*(n-1)-mod(i,2),2),n)];
                b[mod(div(i*(n-1)-mod(i,2),2),n)] = b[mod(div((i-1)*(n-1)-mod(i-1,2),2),n)] - (-1)^i*cos(sum(x[1:i]));
                b[n-mod(div(i*(n-1)-mod(i,2),2),n)] = b[mod(div(i*(n-1)-mod(i,2),2),n)];
            end
        end
        a[m] = 0; b[m] = 1;
        return(a,b)
    else
        throw(DomainError(n,"n must be a power of 2 >= 4"))
    end
end