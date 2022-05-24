% 'cstrt_messine_ngon' provides the vertices coordinates (a,b) of a convex
% small n-gon for n = 2^s and s >= 2
% The 4-gon has maximal width among all 4-gons with fixed perimeter
function [a,b] = cstrt_hansen_ngon(n)
if mod(log2(n),1)==0 && n>=4
    % initialization
    syms u
    v = (pi/2-u)/(n/2-1);
    W = 2*cos(2*u-v)*cos(u)*cos(v-u)/(n*sin(2*u)+(n-2)*sin(2*v-2*u));
    f = diff(W);
    u0 = double(vpasolve(f==0,u,[pi/(2*n-2), pi/n]));
    v0 = (pi/2-u0)/(n/2-1);
    d0 = 2*cos(2*u0-v0)*cos(u0)/(n*sin(2*u0)+(n-2)*sin(2*v0-2*u0));
    D0 = d0*cos(v0-u0)/cos(u0);
    z0 = d0*cos(v0-u0)/cos(2*u0-v0);
    r = ones(n/2+1,1);
    for i = 1:n/4
        r(2*i-1) = d0/D0;
    end
    r(end) = z0/D0;
    x = zeros(n/2,1);
    x(1) = u0;
    for i = 2:n/2
        x(i) = v0;
    end
    % construction
    a = zeros(n-1,1); b = a;
    a(n/2-1) = r(1)*sin(x(1)); a(n/2+1) = -a(n/2-1);
    b(n/2-1) = r(1)*cos(x(1)); b(n/2+1) = b(n/2-1);
    for i = 2:n/2-1
        a(mod(i*(n-1)/2-mod(i,2)/2,n)) = a(mod((i-1)*(n-1)/2-mod(i-1,2)/2,n)) - (-1)^i*r(i)*sin(sum(x(1:i)));
        a(n-mod(i*(n-1)/2-mod(i,2)/2,n)) = -a(mod(i*(n-1)/2-mod(i,2)/2,n));
        b(mod(i*(n-1)/2-mod(i,2)/2,n)) = b(mod((i-1)*(n-1)/2-mod(i-1,2)/2,n)) - (-1)^i*r(i)*cos(sum(x(1:i)));
        b(n-mod(i*(n-1)/2-mod(i,2)/2,n)) = b(mod(i*(n-1)/2-mod(i,2)/2,n));
    end
    a(n/2) = 0; b(n/2) = r(end);
end
end