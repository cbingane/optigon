# 'cstrt_foster_ngon' provides the vertices coordinates (a,b) of a small
# n-gon for even n >= 10
# The 10-gon has maximal area among all small 10-gons
# Please, cite...
#   C. Bingane and M. J. Mossinghoff. Small polygons with large area.
#   arXiv preprint arXiv:2204.04547, 2022.
#   C. Bingane. Tight bounds on the maximal area of small polygons:
#   Improved Mossinghoff polygons. Discrete & Computational Geometry, 2022.
#   M. J. Mossinghoff. Isodiametric problems for polygons. Discrete &
#   Computational Geometry, 36(2): 363-379, 2006.
function cstrt_foster_ngon(n)
    if mod(n,2)==0 && n>=10
        m = div(n,2);
        Pi = big(pi);
        x = Vector{BigFloat}(undef,m);
        # initialization
        v3(u) = (Pi/2-u[1]-2*u[2])/(m-3);
        d3(u) = asin(sin(u[1])-sin(u[1]+u[2]+u[3])+sin(u[1]+2*u[2])+sin(u[1]+2*u[2]+3*v3(u)/2)/(2*cos(v3(u)/2)))-u[1]-2*u[2]-v3(u);
        a30(u) = u[1]; a31(u) = u[2]+u[3]; a32(u) = u[2]-u[3]; a33(u) = v3(u)+d3(u); a34(u) = v3(u)-d3(u);
        p30(u) = a30(u); p31(u) = p30(u)+a31(u); p32(u) = p31(u)+a32(u); p33(u) = p32(u)+a33(u); p34(u) = p33(u)+a34(u);
        x3(u) = sin(p30(u))-sin(p31(u))+sin(p32(u))-sin(p33(u));
        y3(u) = cos(p30(u))-cos(p31(u))+cos(p32(u))-cos(p33(u));
        A31(u) = sin(a30(u)); A32(u) = sin(a31(u)+a32(u))-sin(a31(u));
        A33(u) = (sin(a32(u)+a33(u))-sin(a32(u)))-(sin(a31(u)+a32(u)+a33(u))-sin(a31(u)+a32(u)));
        A34(u) = (sin(a33(u)+a34(u))-sin(a33(u)))-(sin(a32(u)+a33(u)+a34(u))-sin(a32(u)+a33(u)))+(sin(a31(u)+a32(u)+a33(u)+a34(u))-sin(a31(u)+a32(u)+a33(u)));
        A35(u) = (m-5)*(sin(v3(u))-tan(v3(u)/2))-(x3(u)*sin(p34(u))+y3(u)*cos(p34(u))+1/2)*tan(v3(u)/2);
        A3(u) = A31(u)+A32(u)+A33(u)+A34(u)+A35(u);
        F3(u) = -A3(u);
        res = Optim.optimize(F3,[Pi/(2*n-2), Pi/(n-1), 0],LBFGS());
        u0, w10, d10 = Optim.minimizer(res);
        v0 = (Pi/2-u0-2*w10)/(m-3);
        d20 = asin(sin(u0)-sin(u0+w10+d10)+sin(u0+2*w10)+sin(u0+2*w10+3*v0/2)/(2*cos(v0/2)))-u0-2*w10-v0;
        x[1] = u0;
        x[2] = w10+d10; x[3] = w10-d10;
        x[4] = v0+d20; x[5] = v0-d20;
        if n>=12
            x[6:end] .= v0;
        end
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
        throw(DomainError(n,"n must be an even integer >= 10"))
    end
end