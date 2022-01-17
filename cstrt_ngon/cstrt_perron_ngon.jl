# 'cstrt_perron_ngon' provides the vertices coordinates (a,b) of a convex
# equilateral small n-gon for n = 2^s and s >= 4
# Please, cite...
#   C. Bingane and C. Audet. Tight bounds on the maximal perimeter of
#   convex equilateral small polygons. Technical Report G-2021-31, Les
#   cahiers du GERAD, 2021.
function cstrt_perron_ngon(n)
if mod(log2(n),1)==0 && n>=16
    # initialization
    x(t) = sin(t)-(2*cos(t)-1)*cos((n/4-2)*t)*sin(n*t/4)/cos(2*t);
    F(t) = (2*x(t)+sin((n/2-1)*t))^2 + cos((n/2-1)*t)^2 - 4*sin(t/2)^2;
    z = find_zero(F,(2*asin(sin(pi/n)/2),pi/n));
    # construction
    a = zeros(n-1); b = zeros(n-1);
    a[Int(n/2-1)] = sin(z); a[Int(n/2+1)] = -a[Int(n/2-1)];
    b[Int(n/2-1)] = cos(z); b[Int(n/2+1)] = b[Int(n/2-1)];
    a[n-1] = a[Int(n/2-1)]-sin(2*z); a[1] = -a[n-1];
    b[n-1] = b[Int(n/2-1)]-cos(2*z); b[1] = b[n-1];
    a[n-2] = a[Int(n/2-1)]-sin(3*z); a[2] = -a[n-2];
    b[n-2] = b[Int(n/2-1)]-cos(3*z); b[2] = b[n-2];
    for i = 1:Int(n/8-1)
        a[Int(mod(n/2-1+i*(3*n/2-2),n))] = sin(z)-(2*cos(z)-1)*sin(i*(2*z-pi/2))*sin((i+1)*(2*z-pi/2))/cos(2*z);
        a[Int(n-mod(n/2-1+i*(3*n/2-2),n))] = -a[Int(mod(n/2-1+i*(3*n/2-2),n))];
        b[Int(mod(n/2-1+i*(3*n/2-2),n))] = cos(z)-(2*cos(z)-1)*sin(i*(2*z-pi/2))*cos((i+1)*(2*z-pi/2))/cos(2*z);
        b[Int(n-mod(n/2-1+i*(3*n/2-2),n))] = b[Int(mod(n/2-1+i*(3*n/2-2),n))];
        a[Int(mod(n-1+i*(3*n/2-2),n))] = a[Int(mod(n/2-1+i*(3*n/2-2),n))]-(-1)^i*sin((4*i+2)*z);
        a[Int(n-mod(n-1+i*(3*n/2-2),n))] = -a[Int(mod(n-1+i*(3*n/2-2),n))];
        b[Int(mod(n-1+i*(3*n/2-2),n))] = b[Int(mod(n/2-1+i*(3*n/2-2),n))]-(-1)^i*cos((4*i+2)*z);
        b[Int(n-mod(n-1+i*(3*n/2-2),n))] = b[Int(mod(n-1+i*(3*n/2-2),n))];
        a[Int(mod(n-2+i*(3*n/2-2),n))] = a[Int(mod(n/2-1+i*(3*n/2-2),n))]-(-1)^i*sin((4*i+3)*z);
        a[Int(n-mod(n-2+i*(3*n/2-2),n))] = -a[Int(mod(n-2+i*(3*n/2-2),n))];
        b[Int(mod(n-2+i*(3*n/2-2),n))] = b[Int(mod(n/2-1+i*(3*n/2-2),n))]-(-1)^i*cos((4*i+3)*z);
        b[Int(n-mod(n-2+i*(3*n/2-2),n))] = b[Int(mod(n-2+i*(3*n/2-2),n))];
        a[Int(mod(n+i*(3*n/2-2),n))] = a[Int(mod(n/2-1+i*(3*n/2-2),n))]-(-1)^i*sin((4*i+1)*z);
        a[Int(n-mod(n+i*(3*n/2-2),n))] = -a[Int(mod(n+i*(3*n/2-2),n))];
        b[Int(mod(n+i*(3*n/2-2),n))] = b[Int(mod(n/2-1+i*(3*n/2-2),n))]-(-1)^i*cos((4*i+1)*z);
        b[Int(n-mod(n+i*(3*n/2-2),n))] = b[Int(mod(n+i*(3*n/2-2),n))];
    end
    a[Int(n/2)] = 0; b[Int(n/2)] = 1;
    return(a,b)
end
end