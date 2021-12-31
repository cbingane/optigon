# 'cstrt_hansen_ngon' provides the vertices coordinates (a,b) of a convex
# equilateral small n-gon for n = 2^s and s >= 3
function cstrt_hansen_ngon(n)
if mod(log2(n),1)==0 && n>=8
    # initialization
    F(t) = (cos(t)^2+sin((n/2-2)*t)^2-2*cos(t)*sin((n/2-2)*t)^2)/(2*cos(t)^2*(1+sin((n/2-2)*t)))-4*sin(t/2)^2;
    z = find_zero(F,(2*asin(sin(pi/n)/2),pi/n));
    # construction
    a = zeros(n-1); b = zeros(n-1);
    for i = 1:Int(n/4-1)
        a[Int(mod(i*(n/2-1),n))] = (-1)^(i-1)*sin(2*i*z)/2/cos(z);
        a[Int(n-mod(i*(n/2-1),n))] = -a[Int(mod(i*(n/2-1),n))];
        b[Int(mod(i*(n/2-1),n))] = (1+(-1)^(i-1)*cos(2*i*z))/2/cos(z);
        b[Int(n-mod(i*(n/2-1),n))] = b[Int(mod(i*(n/2-1),n))];            
    end
    a[Int(mod((n/4)*(n/2-1),n)+1)] = a[Int(mod((n/4-1)*(n/2-1),n))]+(-1)^(n/4-1)*sin((n/2-2)*z);
    a[Int(n-mod((n/4)*(n/2-1),n)-1)] = -a[Int(mod((n/4)*(n/2-1),n)+1)];
    b[Int(mod((n/4)*(n/2-1),n)+1)] = b[Int(mod((n/4-1)*(n/2-1),n))]+(-1)^(n/4-1)*cos((n/2-2)*z);
    b[Int(n-mod((n/4)*(n/2-1),n)-1)] = b[Int(mod((n/4)*(n/2-1),n)+1)];
    a[Int(n/4)] = 1/2; a[Int(3*n/4)] = -a[Int(n/4)];
    b[Int(n/4)] = (2*cos((n/4-1)*z)-cos(z)*(cos((n/4-1)*z)-sin((n/4-1)*z)))/
    (2*cos(z)*(cos((n/4-1)*z)+sin((n/4-1)*z)));
    b[Int(3*n/4)] = b[Int(n/4)];
    if n>8
        for i = 1:Int(n/4-2)
            a[Int(mod(i*(n/2-1)+n/2,n))] = a[Int(mod(i*(n/2-1),n))] + (-1)^i*sin(2*i*z);
            a[Int(n-mod(i*(n/2-1)+n/2,n))] = -a[Int(mod(i*(n/2-1)+n/2,n))];
            b[Int(mod(i*(n/2-1)+n/2,n))] = b[Int(mod(i*(n/2-1),n))] + (-1)^i*cos(2*i*z);
            b[Int(n-mod(i*(n/2-1)+n/2,n))] = b[Int(mod(i*(n/2-1)+n/2,n))];
        end
    end
    a[Int(n/2)] = 0; b[Int(n/2)] = 1;
    return(a,b)
end
end