# 'cstrt_xiong_ngon' provides the vertices coordinates (a,b) of a small
# n-gon for n = 2m and m >= 4
# The 8-gon has maximal area among all small 8-gons
# Please, cite...
#   C. Bingane and M. J. Mossinghoff. Small polygons with large area.
#   arXiv preprint arXiv:2204.04547, 2022.
#   C. Bingane. Tight bounds on the maximal area of small polygons:
#   Improved Mossinghoff polygons. Discrete & Computational Geometry, 2022.
#   M. J. Mossinghoff. Isodiametric problems for polygons. Discrete &
#   Computational Geometry, 36(2): 363-379, 2006.
function cstrt_xiong_ngon(n)
    if mod(n,2)==0 && n>=8
        m = div(n,2);
        Pi = big(pi);
        x = Vector{BigFloat}(undef,m);
        # initialization
        v2(u) = (Pi/2-u[1]-2*u[2])/(m-3);
        d2(u) = asin(sin(u[1])+sin(u[1]+2*u[2]-v2(u)/2)/(2*cos(v2(u)/2)))-u[1]-u[2];
        A2(u) = sin(u[1])+sin(2*u[2])-sin(u[2]+d2(u))+
        (m-3)*(sin(v2(u))-tan(v2(u)/2))+(cos(u[2]-d2(u))-cos(2*u[2])-1/2)*tan(v2(u)/2);
        F2(u) = -A2(u);
        res = Optim.optimize(F2,[Pi/(2*n-2), Pi/(n-1)],LBFGS());
        u0, w0 = Optim.minimizer(res);
        v0 = (Pi/2-u0-2*w0)/(m-3);
        d0 = asin(sin(u0)+sin(u0+2*w0-v0/2)/(2*cos(v0/2)))-u0-w0;
        x[1] = u0;
        x[2] = w0+d0;
        x[3] = w0-d0;
        x[4:end] .= v0;
        # construction
        a = Vector{BigFloat}(undef,n-1); b = Vector{BigFloat}(undef,n-1);
        a[m-1] = sin(x[1]); a[m+1] = -a[m-1];
        b[m-1] = cos(x[1]); b[m+1] = b[m-1];
        for i = 2:m-1
            a[mod(div(i*(n-1)-mod(i,2),2),n)] = a[mod(div((i-1)*(n-1)-mod(i-1,2),2),n)]-(-1)^i*sin(sum(x[1:i]));
            a[n-mod(div(i*(n-1)-mod(i,2),2),n)] = -a[mod(div(i*(n-1)-mod(i,2),2),n)];
            b[mod(div(i*(n-1)-mod(i,2),2),n)] = b[mod(div((i-1)*(n-1)-mod(i-1,2),2),n)]-(-1)^i*cos(sum(x[1:i]));
            b[n-mod(div(i*(n-1)-mod(i,2),2),n)] = b[mod(div(i*(n-1)-mod(i,2),2),n)];
        end
        a[m] = 0; b[m] = 1;
        return(a,b)
    else
        throw(DomainError(n,"n must be an even integer >= 8"))
    end
end