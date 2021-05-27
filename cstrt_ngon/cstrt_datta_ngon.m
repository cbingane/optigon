function [a,b] = cstrt_datta_ngon(n)
if mod(log(n),log(2))==0 && n>=16
    % initialization
    d = atan(tan(2*pi/n)*tan(pi/n))-asin(sin(2*pi/n)^2/sqrt(2+2*cos(2*pi/n)*cos(4*pi/n)));
    x = zeros(3*n/8,1);
    for i = 1:3*n/8
        if mod(i,3) == 2
            x(i) = 2*(pi/n + (-1)^(i-1)*d);
        else
            x(i) = pi/n + (-1)^(i-1)*d;
        end
    end
    p = zeros(3*n/8-1,1); q = p;
    p(1) = sin(x(1)); q(1) = cos(x(1));
    for i = 2:3*n/8-1
        p(i) = p(i-1) - (-1)^i*sin(sum(x(1:i)));
        q(i) = q(i-1) - (-1)^i*cos(sum(x(1:i)));
    end
    % construction
    a = zeros(n-1,1); b = a;
    a(n/2-1) = p(1); a(n/2+1) = -a(n/2-1);
    b(n/2-1) = q(1); b(n/2+1) = b(n/2-1);
    a(n-1) = a(n/2-1)-sin(x(1)+x(2)/2); a(1) = -a(n-1);
    b(n-1) = b(n/2-1)-cos(x(1)+x(2)/2); b(1) = b(n-1);
    a(n-2) = p(2); a(2) = -a(n-2);
    b(n-2) = q(2); b(2) = b(n-2);
    for i = 1:n/8-1
        a(mod(n/2-1+i*(3*n/2-2),n)) = p(3*i+1);
        a(n-mod(n/2-1+i*(3*n/2-2),n)) = -a(mod(n/2-1+i*(3*n/2-2),n));
        b(mod(n/2-1+i*(3*n/2-2),n)) = q(3*i+1);
        b(n-mod(n/2-1+i*(3*n/2-2),n)) = b(mod(n/2-1+i*(3*n/2-2),n));
        a(mod(n-1+i*(3*n/2-2),n)) = a(mod(n/2-1+i*(3*n/2-2),n))-(-1)^i*sin(sum(x(1:3*i+1))+x(3*i+2)/2);
        a(n-mod(n-1+i*(3*n/2-2),n)) = -a(mod(n-1+i*(3*n/2-2),n));
        b(mod(n-1+i*(3*n/2-2),n)) = b(mod(n/2-1+i*(3*n/2-2),n))-(-1)^i*cos(sum(x(1:3*i+1))+x(3*i+2)/2);
        b(n-mod(n-1+i*(3*n/2-2),n)) = b(mod(n-1+i*(3*n/2-2),n));
        a(mod(n-2+i*(3*n/2-2),n)) = p(3*i+2);
        a(n-mod(n-2+i*(3*n/2-2),n)) = -a(mod(n-2+i*(3*n/2-2),n));
        b(mod(n-2+i*(3*n/2-2),n)) = q(3*i+2);
        b(n-mod(n-2+i*(3*n/2-2),n)) = b(mod(n-2+i*(3*n/2-2),n));
        a(mod(n+i*(3*n/2-2),n)) = p(3*i);
        a(n-mod(n+i*(3*n/2-2),n)) = -a(mod(n+i*(3*n/2-2),n));
        b(mod(n+i*(3*n/2-2),n)) = q(3*i);
        b(n-mod(n+i*(3*n/2-2),n)) = b(mod(n+i*(3*n/2-2),n));
    end
    a(n/2) = 0; b(n/2) = 1;
end
end