% 'cstrt_szabo_ngon' provides the vertices coordinates (a,b) of a small
% n-gon for n = 2m and m >= 6
% The 12-gon has maximal area among all small 12-gons
% Please, cite...
%   C. Bingane and M. J. Mossinghoff. Small polygons with large area.
%   arXiv preprint arXiv:2204.04547, 2022.
%   C. Bingane. Tight bounds on the maximal area of small polygons:
%   Improved Mossinghoff polygons. Discrete & Computational Geometry, 2022.
%   M. J. Mossinghoff. Isodiametric problems for polygons. Discrete &
%   Computational Geometry, 36(2): 363-379, 2006.
function [a,b] = cstrt_szabo_ngon(n)
if mod(n,2) == 0 && n>=12
    % initialization
    syms u w1 w2 d1
    v = (pi/2-u-2*w1-2*w2)/(n/2-5);
    d2 = asin(sin(u)-sin(u+w1+d1)+sin(u+2*w1)+sin(u+2*w1+2*w2-v/2)/(2*cos(v/2)))-u-2*w1-w2;
    a0 = u; a1 = w1+d1; a2 = w1-d1; a3 = w2+d2; a4 = w2-d2;
    p0 = a0; p1 = p0+a1; p2 = p1+a2; p3 = p2+a3; p4 = p3+a4;
    x4 = sin(p0)-sin(p1)+sin(p2)-sin(p3);
    y4 = cos(p0)-cos(p1)+cos(p2)-cos(p3);
    A1 = sin(a0); A2 = sin(a1+a2)-sin(a1);
    A3 = (sin(a2+a3)-sin(a2))-(sin(a1+a2+a3)-sin(a1+a2));
    A4 = (sin(a3+a4)-sin(a3))-(sin(a2+a3+a4)-sin(a2+a3))+(sin(a1+a2+a3+a4)-sin(a1+a2+a3));
    A5 = (n/2-5)*(sin(v)-tan(v/2))-(x4*sin(p4)+y4*cos(p4)+1/2)*tan(v/2);
    A = A1+A2+A3+A4+A5;
    fu = diff(A,u); fw1 = diff(A,w1); fd1 = diff(A,d1); fw2 = diff(A,w2);
    sol = vpasolve([fu==0;fw1==0;fd1==0;fw2==0],[u;w1;d1;w2],[pi/(2*n-2), pi/n; pi/n, 2*pi/n; 0, pi/n; pi/n, 2*pi/n]);
    u0 = double(sol.u); w10 = double(sol.w1); d10 = double(sol.d1); w20 = double(sol.w2);
    v0 = (pi/2-u0-2*w10-2*w20)/(n/2-5);
    d20 = asin(sin(u0)-sin(u0+w10+d10)+sin(u0+2*w10)+sin(u0+2*w10+2*w20-v0/2)/(2*cos(v0/2)))-u0-2*w10-w20;
    x = zeros(n/2,1);
    x(1) = u0;
    x(2) = w10+d10; x(3) = w10-d10;
    x(4) = w20+d20; x(5) = w20-d20;
    x(6:end) = v0;
    % construction
    a = zeros(n-1,1); b = a;
    a(n/2-1) = sin(x(1)); a(n/2+1) = -a(n/2-1);
    b(n/2-1) = cos(x(1)); b(n/2+1) = b(n/2-1);
    for i = 2:n/2-1
        a(mod(i*(n-1)/2-mod(i,2)/2,n)) = a(mod((i-1)*(n-1)/2-mod(i-1,2)/2,n)) - (-1)^i*sin(sum(x(1:i)));
        a(n-mod(i*(n-1)/2-mod(i,2)/2,n)) = -a(mod(i*(n-1)/2-mod(i,2)/2,n));
        b(mod(i*(n-1)/2-mod(i,2)/2,n)) = b(mod((i-1)*(n-1)/2-mod(i-1,2)/2,n)) - (-1)^i*cos(sum(x(1:i)));
        b(n-mod(i*(n-1)/2-mod(i,2)/2,n)) = b(mod(i*(n-1)/2-mod(i,2)/2,n));
    end
    a(n/2) = 0; b(n/2) = 1;
end
end