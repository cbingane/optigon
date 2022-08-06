# 'cstrt_fodor_ngon' provides the vertices coordinates (a,b) of a convex
# equilateral small n-gon for n = 2^s and s >= 3
# The 8-gon has maximal width among all equilateral small 8-gons
# Please, cite...
#   C. Bingane and C. Audet. The equilateral small octagon of maximal
#   width. Mathematics of Computation, 2022.
function cstrt_fodor_ngon(n)
    if mod(log2(n),1)==0 && n>=8
        # initialization
        v(u) = (pi/2-(n/4+1)*u)/(n/4-1);
        d(u) = (cos(v(u))-sin(u)+sin(u+v(u)))/(cos(u)+sin(v(u)));
        F(u) = 1+d(u)^2-2*d(u)*cos(u) - (2-2*cos(v(u)));
        u0 = find_zero(F,(pi/n, 2*pi/(n+4)));
        v0 = (pi/2-(n/4+1)*u0)/(n/4-1);
        d0 = (cos(v0)-sin(u0)+sin(u0+v0))/(cos(u0)+sin(v0));
        r = ones(Int(n/4+1));
        for i = 1:Int(n/8)
            r[2*i] = d0;
        end
        x = zeros(Int(n/4+1));
        x[1] = u0;
        x[2:end-1] .= u0+v0;
        x[end] = u0;
        p = zeros(Int(n/4)); q = zeros(Int(n/4));
        p[1] = r[1]*sin(x[1]); q[1] = r[1]*cos(x[1]);
        for i = 2:Int(n/4)
            p[i] = p[i-1] - (-1)^i*r[i]*sin(sum(x[1:i]));
            q[i] = q[i-1] - (-1)^i*r[i]*cos(sum(x[1:i]));
        end
        # construction
        a = zeros(n-1); b = zeros(n-1);
        a[Int(n/2-1)] = p[1]; a[Int(n/2+1)] = -a[Int(n/2-1)];
        b[Int(n/2-1)] = q[1]; b[Int(n/2+1)] = b[Int(n/2-1)];
        if n>=16
            for i = 2:Int(n/4-1)
                a[Int(mod(i*(n/2-1),n))] = p[i];
                a[Int(n-mod(i*(n/2-1),n))] = -a[Int(mod(i*(n/2-1),n))];
                b[Int(mod(i*(n/2-1),n))] = q[i];
                b[Int(n-mod(i*(n/2-1),n))] = b[Int(mod(i*(n/2-1),n))];
                a[Int(mod((i-1)*(n/2-1)+n/2,n))] = a[Int(mod((i-1)*(n/2-1),n))] - (-1)^i*sin(sum(x[1:i-1])+v0);
                a[Int(n-mod((i-1)*(n/2-1)+n/2,n))] = -a[Int(mod((i-1)*(n/2-1)+n/2,n))];
                b[Int(mod((i-1)*(n/2-1)+n/2,n))] = b[Int(mod((i-1)*(n/2-1),n))] - (-1)^i*cos(sum(x[1:i-1])+v0);
                b[Int(n-mod((i-1)*(n/2-1)+n/2,n))] = b[Int(mod((i-1)*(n/2-1)+n/2,n))];
            end
        end
        a[Int(3*n/4)] = p[Int(n/4)]; a[Int(n/4)] = -a[Int(3*n/4)];
        b[Int(3*n/4)] = q[Int(n/4)]; b[Int(n/4)] = b[Int(3*n/4)];
        a[Int(3*n/4+1)] = a[Int(n/4+1)]-sin(sum(x[1:Int(n/4-1)])+v0); a[Int(n/4-1)] = -a[Int(3*n/4+1)];
        b[Int(3*n/4+1)] = b[Int(n/4+1)]-cos(sum(x[1:Int(n/4-1)])+v0); b[Int(n/4-1)] = b[Int(3*n/4+1)];
        a[Int(n/2)] = 0; b[Int(n/2)] = d0;
        return(a,b)
    else
        throw(DomainError(n,"n must be a power of 2 >= 8"))
    end
end