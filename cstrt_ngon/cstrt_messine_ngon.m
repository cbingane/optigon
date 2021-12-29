function [a,b] = cstrt_messine_ngon(n)
if mod(log2(n),1)==0 && n>=4
    % initialization
    d = pi/4 - asin(cos(pi/n)/sqrt(2));
    x = zeros(n/2,1);
    for i = 1:n/2
        x(i) = pi/n + (-1)^i*d;
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