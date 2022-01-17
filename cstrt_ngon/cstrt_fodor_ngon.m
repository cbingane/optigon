% 'cstrt_fodor_ngon' provides the vertices coordinates (a,b) of a convex
% equilateral small n-gon for n = 2^s and s >= 3
% The 8-gon has maximal width among all equilateral small 8-gons
% Please, cite...
%   C. Bingane and C. Audet. The equilateral small octagon of maximal
%   width. Mathematics of Computation, 2022.
function [a,b] = cstrt_fodor_ngon(n)
if mod(log2(n),1)==0 && n>=8
    % initialization
    syms u
    v = (pi/2-(n/4+1)*u)/(n/4-1);
    d = (cos(v)-sin(u)+sin(u+v))/(cos(u)+sin(v));
    eqn = 1+d^2-2*d*cos(u) == 2-2*cos(v);
    u0 = double(vpasolve(eqn,u,[pi/n, 2*pi/(n+4)]));
    v0 = (pi/2-(n/4+1)*u0)/(n/4-1);
    d0 = (cos(v0)-sin(u0)+sin(u0+v0))/(cos(u0)+sin(v0));
    r = ones(n/4+1,1);
    for i = 1:n/8
        r(2*i) = d0;
    end
    x = zeros(n/4+1,1);
    x(1) = u0;
    x(2:end-1) = u0+v0;
    x(end) = u0;
    p = zeros(n/4,1); q = p;
    p(1) = r(1)*sin(x(1)); q(1) = r(1)*cos(x(1));
    for i = 2:n/4
        p(i) = p(i-1) - (-1)^i*r(i)*sin(sum(x(1:i)));
        q(i) = q(i-1) - (-1)^i*r(i)*cos(sum(x(1:i)));
    end
    % construction
    a = zeros(n-1,1); b = a;
    a(n/2-1) = p(1); a(n/2+1) = -a(n/2-1);
    b(n/2-1) = q(1); b(n/2+1) = b(n/2-1);
    for i = 2:n/4-1
        a(mod(i*(n/2-1),n)) = p(i);
        a(n-mod(i*(n/2-1),n)) = -a(mod(i*(n/2-1),n));
        b(mod(i*(n/2-1),n)) = q(i);
        b(n-mod(i*(n/2-1),n)) = b(mod(i*(n/2-1),n));
        a(mod((i-1)*(n/2-1)+n/2,n)) = a(mod((i-1)*(n/2-1),n)) - (-1)^i*sin(sum(x(1:i-1))+v0);
        a(n-mod((i-1)*(n/2-1)+n/2,n)) = -a(mod((i-1)*(n/2-1)+n/2,n));
        b(mod((i-1)*(n/2-1)+n/2,n)) = b(mod((i-1)*(n/2-1),n)) - (-1)^i*cos(sum(x(1:i-1))+v0);
        b(n-mod((i-1)*(n/2-1)+n/2,n)) = b(mod((i-1)*(n/2-1)+n/2,n));
    end
    a(3*n/4) = p(n/4); a(n/4) = -a(3*n/4);
    b(3*n/4) = q(n/4); b(n/4) = b(3*n/4);
    a(3*n/4+1) = a(n/4+1) - sin(sum(x(1:n/4-1))+v0); a(n/4-1) = -a(3*n/4+1);
    b(3*n/4+1) = b(n/4+1) - cos(sum(x(1:n/4-1))+v0); b(n/4-1) = b(3*n/4+1);
    a(n/2) = 0; b(n/2) = d0;
end
end