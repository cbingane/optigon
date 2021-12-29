function [a,b] = cstrt_bingane_ngon(n)
if mod(log2(n),1)==0 && n>=8
    % initialization
    d = pi/n - asin(sin(2*pi/n)/2);
    % construction
    a = zeros(n-1,1); b = a;
    for i = 1:n/4
        a(mod(i*(n/2-1),n)) = sin(2*i*pi/n)*sin(d + (-1)^(i-1)*pi/n)/sin(2*pi/n);
        a(n-mod(i*(n/2-1),n)) = -a(mod(i*(n/2-1),n));
        b(mod(i*(n/2-1),n)) = (sin(pi/n - d) + cos(2*i*pi/n)*sin(d + (-1)^(i-1)*pi/n))/sin(2*pi/n);
        b(n-mod(i*(n/2-1),n)) = b(mod(i*(n/2-1),n));
    end
    for i = 1:n/4-1
        a(mod(i*(n/2-1)+n/2,n)) = a(mod(i*(n/2-1),n)) + (-1)^i*sin(2*i*pi/n);
        a(n-mod(i*(n/2-1)+n/2,n)) = -a(mod(i*(n/2-1)+n/2,n));
        b(mod(i*(n/2-1)+n/2,n)) = b(mod(i*(n/2-1),n)) + (-1)^i*cos(2*i*pi/n);
        b(n-mod(i*(n/2-1)+n/2,n)) = b(mod(i*(n/2-1)+n/2,n));
    end
    a(n/2) = 0; b(n/2) = 1;
end
end