% 'cstrt_mossinghoff_equi_ngon' provides the vertices coordinates (a,b) of
% the Mossinghoff equilateral n-gon for the maximal perimeter problem
% Please, cite...
%   M. J. Mossinghoff. An isodiametric problem for equilateral polygons.
%   Contemporary Mathematics, 457:237-252, 2008.
function [a,b] = cstrt_mossinghoff_equi_ngon(n)
if mod(log(n),log(2))==0 && n>=16
    % initialization
    syms t
    x = sin(t)-sin(2*t)+(2*cos(t)-1)*(sin(2*t)+sin((n/2-2)*t))/(2*cos(2*t));
    eqn = (2*x-sin((n/2-1)*t))^2 + cos((n/2-1)*t)^2 == 4*sin(t/2)^2;
    s = vpasolve(eqn,t,[0,pi/n]);
    z = double(s);
    % construction
    a = zeros(n-1,1); b = a;
    a(n/2-1) = sin(z); a(n/2+1) = -a(n/2-1);
    b(n/2-1) = cos(z); b(n/2+1) = b(n/2-1);
    a(n-1) = a(n/2-1)-sin(2*z); a(1) = -a(n-1);
    b(n-1) = b(n/2-1)-cos(2*z); b(1) = b(n-1);
    a(n/2-2) = a(n-1)+sin(3*z); a(n/2+2) = -a(n/2-2);
    b(n/2-2) = b(n-1)+cos(3*z); b(n/2+2) = b(n/2-2);
    for i = 1:n/8-1
       a(mod(n/2-3+(i-1)*(3*n/2-2),n)) = sin(z)-sin(2*z)+(2*cos(z)-1)*sin(i*(2*z-pi/2))*sin((i+1)*(2*z-pi/2))/cos(2*z);
       a(n-mod(n/2-3+(i-1)*(3*n/2-2),n)) = -a(mod(n/2-3+(i-1)*(3*n/2-2),n));
       b(mod(n/2-3+(i-1)*(3*n/2-2),n)) = cos(z)-cos(2*z)+(2*cos(z)-1)*sin(i*(2*z-pi/2))*cos((i+1)*(2*z-pi/2))/cos(2*z);
       b(n-mod(n/2-3+(i-1)*(3*n/2-2),n)) = b(mod(n/2-3+(i-1)*(3*n/2-2),n));
       a(mod(n-3+(i-1)*(3*n/2-2),n)) = a(mod(n/2-3+(i-1)*(3*n/2-2),n))+(-1)^i*sin((4*i+2)*z);
       a(n-mod(n-3+(i-1)*(3*n/2-2),n)) = -a(mod(n-3+(i-1)*(3*n/2-2),n));
       b(mod(n-3+(i-1)*(3*n/2-2),n)) = b(mod(n/2-3+(i-1)*(3*n/2-2),n))+(-1)^i*cos((4*i+2)*z);
       b(n-mod(n-3+(i-1)*(3*n/2-2),n)) = b(mod(n-3+(i-1)*(3*n/2-2),n));
       a(mod(n-4+(i-1)*(3*n/2-2),n)) = a(mod(n/2-3+(i-1)*(3*n/2-2),n))+(-1)^i*sin((4*i+3)*z);
       a(n-mod(n-4+(i-1)*(3*n/2-2),n)) = -a(mod(n-4+(i-1)*(3*n/2-2),n));
       b(mod(n-4+(i-1)*(3*n/2-2),n)) = b(mod(n/2-3+(i-1)*(3*n/2-2),n))+(-1)^i*cos((4*i+3)*z);
       b(n-mod(n-4+(i-1)*(3*n/2-2),n)) = b(mod(n-4+(i-1)*(3*n/2-2),n));
       a(mod(n-2+(i-1)*(3*n/2-2),n)) = a(mod(n/2-3+(i-1)*(3*n/2-2),n))+(-1)^i*sin((4*i+1)*z);
       a(n-mod(n-2+(i-1)*(3*n/2-2),n)) = -a(mod(n-2+(i-1)*(3*n/2-2),n));
       b(mod(n-2+(i-1)*(3*n/2-2),n)) = b(mod(n/2-3+(i-1)*(3*n/2-2),n))+(-1)^i*cos((4*i+1)*z);
       b(n-mod(n-2+(i-1)*(3*n/2-2),n)) = b(mod(n-2+(i-1)*(3*n/2-2),n));
    end
    a(n/2) = 0; b(n/2) = 1;
end
end