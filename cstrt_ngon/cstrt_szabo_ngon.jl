# 'cstrt_szabo_ngon' provides the vertices coordinates (a,b) of a small
# n-gon for even n >= 12
# The 12-gon has maximal area among all small 12-gons
# Please, cite...
#   C. Bingane and M. J. Mossinghoff. Small polygons with large area.
#   arXiv preprint arXiv:2204.04547, 2022.
#   C. Bingane. Tight bounds on the maximal area of small polygons:
#   Improved Mossinghoff polygons. Discrete & Computational Geometry, 2022.
#   M. J. Mossinghoff. Isodiametric problems for polygons. Discrete &
#   Computational Geometry, 36(2): 363-379, 2006.
function cstrt_szabo_ngon(n)
    if mod(n,2)==0 && n>=12
        m = div(n,2);
        Pi = big(pi);
        x = Vector{BigFloat}(undef,m);
        # initialization
        v4(u) = (Pi/2-u[1]-2*u[2]-2*u[4])/(m-5);
        d4(u) = asin(sin(u[1])-sin(u[1]+u[2]+u[3])+sin(u[1]+2*u[2])+sin(u[1]+2*u[2]+2*u[4]-v4(u)/2)/(2*cos(v4(u)/2)))-u[1]-2*u[2]-u[4];
        a40(u) = u[1]; a41(u) = u[2]+u[3]; a42(u) = u[2]-u[3]; a43(u) = u[4]+d4(u); a44(u) = u[4]-d4(u);
        p40(u) = a40(u); p41(u) = p40(u)+a41(u); p42(u) = p41(u)+a42(u); p43(u) = p42(u)+a43(u); p44(u) = p43(u)+a44(u);
        x4(u) = sin(p40(u))-sin(p41(u))+sin(p42(u))-sin(p43(u));
        y4(u) = cos(p40(u))-cos(p41(u))+cos(p42(u))-cos(p43(u));
        A41(u) = sin(a40(u)); A42(u) = sin(a41(u)+a42(u))-sin(a41(u));
        A43(u) = (sin(a42(u)+a43(u))-sin(a42(u)))-(sin(a41(u)+a42(u)+a43(u))-sin(a41(u)+a42(u)));
        A44(u) = (sin(a43(u)+a44(u))-sin(a43(u)))-(sin(a42(u)+a43(u)+a44(u))-sin(a42(u)+a43(u)))+(sin(a41(u)+a42(u)+a43(u)+a44(u))-sin(a41(u)+a42(u)+a43(u)));
        A45(u) = (m-5)*(sin(v4(u))-tan(v4(u)/2))-(x4(u)*sin(p44(u))+y4(u)*cos(p44(u))+1/2)*tan(v4(u)/2);
        A4(u) = A41(u)+A42(u)+A43(u)+A44(u)+A45(u);
        F4(u) = -A4(u);
        res = Optim.optimize(F4,[Pi/(2*n-2), Pi/(n-1), 0, Pi/(n-1)],LBFGS());
        u0, w10, d10, w20 = Optim.minimizer(res);
        v0 = (Pi/2-u0-2*w10-2*w20)/(m-5);
        d20 = asin(sin(u0)-sin(u0+w10+d10)+sin(u0+2*w10)+sin(u0+2*w10+2*w20-v0/2)/(2*cos(v0/2)))-u0-2*w10-w20;
        x[1] = u0;
        x[2] = w10+d10; x[3] = w10-d10;
        x[4] = w20+d20; x[5] = w20-d20;
        x[6:end] .= v0;
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
        throw(DomainError(n,"n must be an even integer >= 12"))
    end
end