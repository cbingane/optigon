% 'cstrt_foster_ngon' provides the vertices coordinates (a,b) of a small
% n-gon for n = 2m and m >= 3
function [a,b] = cstrt_foster_ngon(n)
if mod(n,2) == 0 && n>=6
    % initialization
    s = (2*sqrt(114)-7)/22;
    u = s*pi/n;
    v = (pi/2-u)/(n/2-1);
    w = asin(sin(u)+sin(u+3*v/2)/(2*cos(v/2)))-u-v;
    x = zeros(n/2,1);
    x(1) = u;
    x(2) = v+w;
    x(3) = v-w;
    if n >= 8
        x(4:end) = v;
    end
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