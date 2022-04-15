% 'cstrt_xiong_ngon' provides the vertices coordinates (a,b) of a small
% n-gon for n = 2m and m >= 4
% The 8-gon has maximal area among all small 8-gons
% Please, cite...
%   C. Bingane and M. J. Mossinghoff. Small polygons with large area.
%   arXiv preprint arXiv:2204.04547, 2022.
%   C. Bingane. Tight bounds on the maximal area of small polygons:
%   Improved Mossinghoff polygons. Discrete & Computational Geometry, 2022.
%   M. J. Mossinghoff. Isodiametric problems for polygons. Discrete &
%   Computational Geometry, 36(2): 363-379, 2006.
function [a,b] = cstrt_xiong_ngon(n)
if mod(n,2) == 0 && n>=8
    % initialization
    syms u w
    v = (pi/2-u-2*w)/(n/2-3);
    d = asin(sin(u)+sin(u+2*w-v/2)/(2*cos(v/2)))-u-w;
    A = sin(u)+sin(2*w)-sin(w+d)+(n/2-3)*(sin(v)-tan(v/2))+(cos(w-d)-cos(2*w)-1/2)*tan(v/2);
    fu = diff(A,u); fw = diff(A,w);
    sol = vpasolve([fu==0;fw==0],[u;w],[pi/(2*n-2), pi/n; pi/n, 2*pi/n]);
    u0 = double(sol.u); w0 = double(sol.w);
    v0 = (pi/2-u0-2*w0)/(n/2-3);
    d0 = asin(sin(u0)+sin(u0+2*w0-v0/2)/(2*cos(v0/2)))-u0-w0;
    x = zeros(n/2,1);
    x(1) = u0;
    x(2) = w0+d0;
    x(3) = w0-d0;
    x(4:end) = v0;
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