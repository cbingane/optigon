% 'cstrt_vincze_ngon' provides the vertices coordinates (a,b) of a convex
% equilateral small n-gon for n = 2^s and s >= 3
% The 8-gon has maximal perimeter among all convex equilateral small 8-gons
% Please, cite...
%   C. Bingane and C. Audet. Tight bounds on the maximal perimeter of
%   convex equilateral small polygons. Archiv der Mathematik, 2022.
%   C. Audet, P. Hansen, F. Messine, and S. Perron. The minimum diameter
%   octagon with unit-length sides: Vincze's wife octagon is suboptimal.
%   Journal of Combinatorial Theory, Series A, 108(1): 63-75, 2004.
function [a,b] = cstrt_vincze_ngon(n)
if mod(log2(n),1)==0 && n>=8
    % initialization
    syms t
    eqn = (cos(t)^2+sin((n/2-2)*t)^2-2*cos(t)*sin((n/2-2)*t)^2)/(2*cos(t)^2*(1+sin((n/2-2)*t))) == 4*sin(t/2)^2;
    z = double(vpasolve(eqn,t,[0,pi/n]));
    % construction
    a = zeros(n-1,1); b = a;
    for i = 1:n/4-1
        a(mod(i*(n/2-1),n)) = (-1)^(i-1)*sin(2*i*z)/2/cos(z); a(n-mod(i*(n/2-1),n)) = -a(mod(i*(n/2-1),n));
        b(mod(i*(n/2-1),n)) = (1+(-1)^(i-1)*cos(2*i*z))/2/cos(z); b(n-mod(i*(n/2-1),n)) = b(mod(i*(n/2-1),n));            
    end
    a(mod((n/4)*(n/2-1),n)+1) = a(mod((n/4-1)*(n/2-1),n))+(-1)^(n/4-1)*sin((n/2-2)*z);
    a(n-mod((n/4)*(n/2-1),n)-1) = -a(mod((n/4)*(n/2-1),n)+1);
    b(mod((n/4)*(n/2-1),n)+1) = b(mod((n/4-1)*(n/2-1),n))+(-1)^(n/4-1)*cos((n/2-2)*z);
    b(n-mod((n/4)*(n/2-1),n)-1) = b(mod((n/4)*(n/2-1),n)+1);
    a(n/4) = 1/2; a(3*n/4) = -a(n/4);
    b(n/4) = (2*cos((n/4-1)*z)-cos(z)*(cos((n/4-1)*z)-sin((n/4-1)*z)))/(2*cos(z)*(cos((n/4-1)*z)+sin((n/4-1)*z)));
    b(3*n/4) = b(n/4);
    if n>8
        for i = 1:n/4-2
            a(mod(i*(n/2-1)+n/2,n)) = a(mod(i*(n/2-1),n)) + (-1)^i*sin(2*i*z);
            a(n-mod(i*(n/2-1)+n/2,n)) = -a(mod(i*(n/2-1)+n/2,n));
            b(mod(i*(n/2-1)+n/2,n)) = b(mod(i*(n/2-1),n)) + (-1)^i*cos(2*i*z);
            b(n-mod(i*(n/2-1)+n/2,n)) = b(mod(i*(n/2-1)+n/2,n));
        end
    end
    a(n/2) = 0; b(n/2) = 1;
end
end