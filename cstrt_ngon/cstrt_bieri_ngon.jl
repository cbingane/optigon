# 'cstrt_bieri_ngon' provides the vertices coordinates (a,b) of a small
# n-gon for n = 2m and m >= 3
# The 6-gon has maximal area among all small 6-gons
# Please, cite...
#   C. Bingane. Tight bounds on the maximal area of small polygons:
#   Improved Mossinghoff polygons. Discrete & Computational Geometry, 2022.
#   M. J. Mossinghoff. Isodiametric problems for polygons. Discrete &
#   Computational Geometry, 36(2): 363-379, 2006.
function cstrt_bieri_ngon(n)
    if mod(n,2)==0 && n>=6
        m = div(n,2);
        Pi = big(pi);
        x = Vector{BigFloat}(undef,m);
        # initialization
        v1(u) = (Pi/2-u)/(m-1);
        d1(u) = asin(sin(u)+sin(u+3*v1(u)/2)/(2*cos(v1(u)/2)))-u-v1(u);
        A1(u) = sin(u)+sin(2*v1(u))-sin(v1(u)+d1(u))+
        (m-3)*(sin(v1(u))-tan(v1(u)/2))+(cos(v1(u)-d1(u))-cos(2*v1(u))-1/2)*tan(v1(u)/2);
        F1(u) = -A1(u);
        res = Optim.optimize(F1,Pi/(2*n-2),Pi/n);
        u0 = Optim.minimizer(res);
        v0 = (Pi/2-u0)/(m-1);
        d0 = asin(sin(u0)+sin(u0+3*v0/2)/(2*cos(v0/2)))-u0-v0;
        x[1] = u0;
        x[2] = v0+d0; x[3] = v0-d0;
        if n>=8
            x[4:end] .= v0;
        end
        # construction
        a = Vector{BigFloat}(undef,n-1); b = Vector{BigFloat}(undef,n-1);
        a[m-1] = sin(x[1]); a[m+1] = -a[m-1];
        b[m-1] = cos(x[1]); b[m+1] = b[m-1];
        if n>=6
            for i = 2:m-1
                a[mod(div(i*(n-1)-mod(i,2),2),n)] = a[mod(div((i-1)*(n-1)-mod(i-1,2),2),n)]-(-1)^i*sin(sum(x[1:i]));
                a[n-mod(div(i*(n-1)-mod(i,2),2),n)] = -a[mod(div(i*(n-1)-mod(i,2),2),n)];
                b[mod(div(i*(n-1)-mod(i,2),2),n)] = b[mod(div((i-1)*(n-1)-mod(i-1,2),2),n)]-(-1)^i*cos(sum(x[1:i]));
                b[n-mod(div(i*(n-1)-mod(i,2),2),n)] = b[mod(div(i*(n-1)-mod(i,2),2),n)];
            end
        end
        a[m] = 0; b[m] = 1;
        return(a,b)
    else
        throw(DomainError(n,"n must be an even integer >= 6"))
    end
end