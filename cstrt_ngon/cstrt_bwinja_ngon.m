% 'cstrt_bwinja_ngon' provides the vertices coordinates (a,b) of a small
% n-gon for n = 2m and m >= 3
function [a,b] = cstrt_bwinja_ngon(n)
if mod(n,2) == 0 && n>=6
    % initialization
    syms u
    v = (pi/2-u)/(n/2-1);
    d = asin(sin(u)+sin(u+3*v/2)/(2*cos(v/2)))-u-v;
    A = sin(u)+sin(2*v)-sin(v+d)+...
        (n/2-3)*(sin(v)-tan(v/2))+(cos(v-d)-cos(2*v)-1/2)*tan(v/2);
    f = diff(A);
    u0 = double(vpasolve(f==0,u,[pi/(2*n-2), pi/n]));
    v0 = (pi/2-u0)/(n/2-1);
    d0 = asin(sin(u0)+sin(u0+3*v0/2)/(2*cos(v0/2)))-u0-v0;
    x = zeros(n/2,1);
    x(1) = u0;
    x(2) = v0+d0;
    x(3) = v0-d0;
    if n >= 8
        x(4:end) = v0;
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