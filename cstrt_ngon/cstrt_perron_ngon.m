% 'cstrt_audet_ngon' provides the vertices coordinates (a,b) of a convex
% equilateral small n-gon for n = 2^s and s >= 4
% Please, cite...
%   C. Bingane and C. Audet. Tight bounds on the maximal perimeter of
%   convex equilateral small polygons. Technical Report G-2021-31, Les
%   cahiers du GERAD, 2021.
function [a,b] = cstrt_perron_ngon(n)
if mod(log2(n),1)==0 && n>=16
    % initialization
    syms t
    x = sin(t)-(2*cos(t)-1)*cos((n/4-2)*t)*sin(n*t/4)/cos(2*t);
    eqn = (2*x+sin((n/2-1)*t))^2 + cos((n/2-1)*t)^2 == 4*sin(t/2)^2;
    z = double(vpasolve(eqn,t,[0,pi/n]));
    % construction
    a = zeros(n-1,1); b = a;
    a(n/2-1) = sin(z); a(n/2+1) = -a(n/2-1);
    b(n/2-1) = cos(z); b(n/2+1) = b(n/2-1);
    a(n-1) = a(n/2-1)-sin(2*z); a(1) = -a(n-1);
    b(n-1) = b(n/2-1)-cos(2*z); b(1) = b(n-1);
    a(n-2) = a(n/2-1)-sin(3*z); a(2) = -a(n-2);
    b(n-2) = b(n/2-1)-cos(3*z); b(2) = b(n-2);
    for i = 1:n/8-1
        a(mod(n/2-1+i*(3*n/2-2),n)) = sin(z)-(2*cos(z)-1)*sin(i*(2*z-pi/2))*sin((i+1)*(2*z-pi/2))/cos(2*z);
        a(n-mod(n/2-1+i*(3*n/2-2),n)) = -a(mod(n/2-1+i*(3*n/2-2),n));
        b(mod(n/2-1+i*(3*n/2-2),n)) = cos(z)-(2*cos(z)-1)*sin(i*(2*z-pi/2))*cos((i+1)*(2*z-pi/2))/cos(2*z);
        b(n-mod(n/2-1+i*(3*n/2-2),n)) = b(mod(n/2-1+i*(3*n/2-2),n));
        a(mod(n-1+i*(3*n/2-2),n)) = a(mod(n/2-1+i*(3*n/2-2),n))-(-1)^i*sin((4*i+2)*z);
        a(n-mod(n-1+i*(3*n/2-2),n)) = -a(mod(n-1+i*(3*n/2-2),n));
        b(mod(n-1+i*(3*n/2-2),n)) = b(mod(n/2-1+i*(3*n/2-2),n))-(-1)^i*cos((4*i+2)*z);
        b(n-mod(n-1+i*(3*n/2-2),n)) = b(mod(n-1+i*(3*n/2-2),n));
        a(mod(n-2+i*(3*n/2-2),n)) = a(mod(n/2-1+i*(3*n/2-2),n))-(-1)^i*sin((4*i+3)*z);
        a(n-mod(n-2+i*(3*n/2-2),n)) = -a(mod(n-2+i*(3*n/2-2),n));
        b(mod(n-2+i*(3*n/2-2),n)) = b(mod(n/2-1+i*(3*n/2-2),n))-(-1)^i*cos((4*i+3)*z);
        b(n-mod(n-2+i*(3*n/2-2),n)) = b(mod(n-2+i*(3*n/2-2),n));
        a(mod(n+i*(3*n/2-2),n)) = a(mod(n/2-1+i*(3*n/2-2),n))-(-1)^i*sin((4*i+1)*z);
        a(n-mod(n+i*(3*n/2-2),n)) = -a(mod(n+i*(3*n/2-2),n));
        b(mod(n+i*(3*n/2-2),n)) = b(mod(n/2-1+i*(3*n/2-2),n))-(-1)^i*cos((4*i+1)*z);
        b(n-mod(n+i*(3*n/2-2),n)) = b(mod(n+i*(3*n/2-2),n));
    end
    a(n/2) = 0; b(n/2) = 1;
end
end