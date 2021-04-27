% 'cstrt_audet_ngon' provides the vertices coordinates (a,b) of a convex
% equilateral small n-gon for n = 2^s and s >= 4
function [a,b] = cstrt_audet_ngon(n)
if mod(log(n),log(2))==0 && n>=16
    % initialization
    syms t x y
    eqn1 = x == 1/2;
    eqn2 = y == 1/2-(1-cos(t))*(cos((n/4-1)*t)-sin((n/4-1)*t)-sin(t))/(cos(t)*(cos(t)+sin(t)));
    eqn3 = (x-(1+cos((n/4-2)*t)-sin((n/4)*t))/(2*cos(t)))^2 + (y-(1+sin((n/4-2)*t)-cos((n/4)*t))/(2*cos(t)))^2 == 4*sin(t/2)^2;
    s = vpasolve([eqn1;eqn2;eqn3],[t;x;y],[0,pi/n;-1/2,1/2;0,1]);
    z = double(s.t);
    % construction
    a = zeros(n-1,1); b = a;
    for i = 1:n/8
        a(mod(i*(n/2-1),n)) = (-1)^(i-1)*sin(2*i*z)/2/cos(z); a(n-mod(i*(n/2-1),n)) = -a(mod(i*(n/2-1),n));
        b(mod(i*(n/2-1),n)) = (1+(-1)^(i-1)*cos(2*i*z))/2/cos(z); b(n-mod(i*(n/2-1),n)) = b(mod(i*(n/2-1),n));            
    end
    for i = 1:n/8
        a(mod(i*(n/2-1)+n/2,n)) = a(mod(i*(n/2-1),n)) + (-1)^i*sin(2*i*z);
        a(n-mod(i*(n/2-1)+n/2,n)) = -a(mod(i*(n/2-1)+n/2,n));
        b(mod(i*(n/2-1)+n/2,n)) = b(mod(i*(n/2-1),n)) + (-1)^i*cos(2*i*z);
        b(n-mod(i*(n/2-1)+n/2,n)) = b(mod(i*(n/2-1)+n/2,n));
    end
    a(n/4) = double(s.x); a(3*n/4) = -a(n/4);
    b(n/4) = double(s.y); b(3*n/4) = b(n/4);
    for i = 1:n/8-1
        a(mod(i*(n/2+1)+n/4,n)) = a(n/4) - (1+(-1)^(i-1)*cos(2*i*z))/2/cos(z);
        a(n-mod(i*(n/2+1)+n/4,n)) = -a(mod(i*(n/2+1)+n/4,n));
        b(mod(i*(n/2+1)+n/4,n)) = b(n/4) - (-1)^(i-1)*sin(2*i*z)/2/cos(z);
        b(n-mod(i*(n/2+1)+n/4,n)) = b(mod(i*(n/2+1)+n/4,n));
    end
    for i = 1:n/8-1
        a(mod(i*(n/2+1)-n/4,n)) = a(mod(i*(n/2+1)+n/4,n)) - (-1)^i*cos(2*i*z);
        a(n-mod(i*(n/2+1)-n/4,n)) = -a(mod(i*(n/2+1)-n/4,n));
        b(mod(i*(n/2+1)-n/4,n)) = b(mod(i*(n/2+1)+n/4,n)) - (-1)^i*sin(2*i*z);
        b(n-mod(i*(n/2+1)-n/4,n)) = b(mod(i*(n/2+1)-n/4,n));
    end
    a(n/2) = 0; b(n/2) = 1;
end
end