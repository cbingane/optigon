% 'cstrt_bezdek_ngon' provides the vertices coordinates (a,b) of a convex
% equilateral small n-gon for n = 2^s and s >= 4
% Please, cite...
%   C. Bingane and C. Audet. The equilateral small octagon of maximal
%   width. Mathematics of Computation, 2022.
function [a,b] = cstrt_bezdek_ngon(n)
if mod(log2(n),1)==0 && n>=16
    % initialization
    syms u
    v = 4*pi/n-3*u;
    d = (2*cos(2*u+v)+1)/(2*cos(u)+cos(3*u+v));
    eqn = 1+d^2-2*d*cos(u) == 2-2*cos(v);
    u0 = double(vpasolve(eqn,u,[pi/n, 4*pi/(3*n)]));
    v0 = 4*pi/n-3*u0;
    d0 = (2*cos(2*u0+v0)+1)/(2*cos(u0)+cos(3*u0+v0));
    r = ones(3*n/8,1);
    for i = 1:3*n/16
        r(2*i) = d0;
    end
    x = zeros(3*n/8,1);
    for i = 1:3*n/8
        if mod(i,3) == 2
            x(i) = u0+v0;
        else
            x(i) = u0;
        end
    end
    p = zeros(3*n/8-1,1); q = p;
    p(1) = r(1)*sin(x(1)); q(1) = r(1)*cos(x(1));
    for i = 2:3*n/8-1
        p(i) = p(i-1) - (-1)^i*r(i)*sin(sum(x(1:i)));
        q(i) = q(i-1) - (-1)^i*r(i)*cos(sum(x(1:i)));
    end
    % construction
    a = zeros(n-1,1); b = a;
    a(n/2-1) = p(1); a(n/2+1) = -a(n/2-1);
    b(n/2-1) = q(1); b(n/2+1) = b(n/2-1);
    a(n-1) = a(n/2-1)-sin(x(1)+v0); a(1) = -a(n-1);
    b(n-1) = b(n/2-1)-cos(x(1)+v0); b(1) = b(n-1);
    a(n-2) = p(2); a(2) = -a(n-2);
    b(n-2) = q(2); b(2) = b(n-2);
    for i = 1:n/8-1
        a(mod(n/2-1+i*(3*n/2-2),n)) = p(3*i+1);
        a(n-mod(n/2-1+i*(3*n/2-2),n)) = -a(mod(n/2-1+i*(3*n/2-2),n));
        b(mod(n/2-1+i*(3*n/2-2),n)) = q(3*i+1);
        b(n-mod(n/2-1+i*(3*n/2-2),n)) = b(mod(n/2-1+i*(3*n/2-2),n));
        a(mod(n-1+i*(3*n/2-2),n)) = a(mod(n/2-1+i*(3*n/2-2),n))-(-1)^i*sin(sum(x(1:3*i+1))+u0);
        a(n-mod(n-1+i*(3*n/2-2),n)) = -a(mod(n-1+i*(3*n/2-2),n));
        b(mod(n-1+i*(3*n/2-2),n)) = b(mod(n/2-1+i*(3*n/2-2),n))-(-1)^i*cos(sum(x(1:3*i+1))+u0);
        b(n-mod(n-1+i*(3*n/2-2),n)) = b(mod(n-1+i*(3*n/2-2),n));
        a(mod(n-2+i*(3*n/2-2),n)) = p(3*i+2);
        a(n-mod(n-2+i*(3*n/2-2),n)) = -a(mod(n-2+i*(3*n/2-2),n));
        b(mod(n-2+i*(3*n/2-2),n)) = q(3*i+2);
        b(n-mod(n-2+i*(3*n/2-2),n)) = b(mod(n-2+i*(3*n/2-2),n));
        a(mod(n+i*(3*n/2-2),n)) = p(3*i);
        a(n-mod(n+i*(3*n/2-2),n)) = -a(mod(n+i*(3*n/2-2),n));
        b(mod(n+i*(3*n/2-2),n)) = q(3*i);
        b(n-mod(n+i*(3*n/2-2),n)) = b(mod(n+i*(3*n/2-2),n));
    end
    a(n/2) = 0; b(n/2) = d0;
end
end