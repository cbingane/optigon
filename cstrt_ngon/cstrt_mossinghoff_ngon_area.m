% 'cstrt_mossinghoff_ngon_area' provides the vertices coordinates (a,b) of
% the Mossinghoff n-gon for the maximal area problem
% Please, cite...
%   M. J. Mossinghoff. Isodiametric problems for polygons. Discrete &
%   Computational Geometry, 36(2): 363-379, 2006.
function [a,b] = cstrt_mossinghoff_ngon_area(n)
if mod(n,2) == 0 && n>=6
    % initialization
    x = zeros(n/2-2,1);
    s = (2*sqrt(114)-7)/22;
    t = (44*(103104*sqrt(114)-998743)+(-1)^(n/2)*75*pi*(347*sqrt(114)-714))/8811220;
    x(1) = s*pi/n + t*pi/n^2;
    if n >= 8
        u = 2*(1-s); v = (2*s-1)/4; w = (s+t-1)/2;
        x(2) = (1+v)*pi/n + (u+w)*pi/n^2;
        if n >= 10
            x(3) = (1-v)*pi/n + (u-w)*pi/n^2;
            if n >= 12
                x(4:end) = pi/n + u*pi/n^2;
            end
        end
    end
    % construction
    a = zeros(n-1,1); b = a;
    a(n/2-1) = sin(x(1)); a(n/2+1) = -a(n/2-1);
    b(n/2-1) = cos(x(1)); b(n/2+1) = b(n/2-1);
    if n >= 8
        for i = 2:n/2-2
            a(mod(i*(n-1)/2-mod(i,2)/2,n)) = a(mod((i-1)*(n-1)/2-mod(i-1,2)/2,n)) - (-1)^i*sin(sum(x(1:i)));
            a(n-mod(i*(n-1)/2-mod(i,2)/2,n)) = -a(mod(i*(n-1)/2-mod(i,2)/2,n));
            b(mod(i*(n-1)/2-mod(i,2)/2,n)) = b(mod((i-1)*(n-1)/2-mod(i-1,2)/2,n)) - (-1)^i*cos(sum(x(1:i)));
            b(n-mod(i*(n-1)/2-mod(i,2)/2,n)) = b(mod(i*(n-1)/2-mod(i,2)/2,n));
        end
    end
    a(mod((n/2-1)*(n-1)/2-mod(n/2-1,2)/2,n)) = -(-1)^(n/2-1)/2;
    a(n-mod((n/2-1)*(n-1)/2-mod(n/2-1,2)/2,n)) = -a(mod((n/2-1)*(n-1)/2-mod(n/2-1,2)/2,n));
    b(mod((n/2-1)*(n-1)/2-mod(n/2-1,2)/2,n)) = b(mod((n/2-2)*(n-1)/2-mod(n/2-2,2)/2,n)) - (-1)^(n/2-1)*sqrt(1 - (a(mod((n/2-1)*(n-1)/2-mod(n/2-1,2)/2,n)) - a(mod((n/2-2)*(n-1)/2-mod(n/2-2,2)/2,n)))^2);
    b(n-mod((n/2-1)*(n-1)/2-mod(n/2-1,2)/2,n)) = b(mod((n/2-1)*(n-1)/2-mod(n/2-1,2)/2,n));
    a(n/2) = 0; b(n/2) = 1;
end
end