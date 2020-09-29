function [a,b] = cstrt_mossinghoff_ngon(n)
if mod(log(n),log(2))==0 && n>=8
    % initialization
    x = pi/n + pi^2/(2*n^2) - pi^2/n^3;
    y = pi/n - pi^2/n^3;
    % construction
    a = zeros(n-1,1); b = a;
    for i = 1:n/4-1
        a(mod(i*(n/2-1),n)) = (sin(x-y) + (-1)^(i-1)*sin(x+(2*i-1)*y))/(2*cos(y));
        a(n-mod(i*(n/2-1),n)) = -a(mod(i*(n/2-1),n));
        b(mod(i*(n/2-1),n)) = (cos(x-y) + (-1)^(i-1)*cos(x+(2*i-1)*y))/(2*cos(y));
        b(n-mod(i*(n/2-1),n)) = b(mod(i*(n/2-1),n));
    end
    if n>8
        for i = 1:n/4-2
            a(mod(i*(n/2-1)+n/2,n)) = a(mod(i*(n/2-1),n)) + (-1)^i*sin(x+(2*i-1)*y);
            a(n-mod(i*(n/2-1)+n/2,n)) = -a(mod(i*(n/2-1)+n/2,n));
            b(mod(i*(n/2-1)+n/2,n)) = b(mod(i*(n/2-1),n)) + (-1)^i*cos(x+(2*i-1)*y);
            b(n-mod(i*(n/2-1)+n/2,n)) = b(mod(i*(n/2-1)+n/2,n));
        end
    end
    a(3*n/4) = -1/2; a(n/4) = -a(3*n/4);
    b(3*n/4) = b(n/4+1) - sqrt(1 - (a(3*n/4)-a(n/4+1))^2); b(n/4) = b(3*n/4);
    if n == 8
        z = asin(abs(a(3*n/4)+1i*b(3*n/4))/2);
    else
        z = asin(abs((a(3*n/4)-a(3*n/4+2))+1i*(b(3*n/4)-b(3*n/4+2)))/2);
    end
    a(3*n/4+1) = a(n/4+1) - sin(x+(n/4-2)*2*y+z); a(n/4-1) = -a(3*n/4+1);
    b(3*n/4+1) = b(n/4+1) - cos(x+(n/4-2)*2*y+z); b(n/4-1) = b(3*n/4+1);
    a(n/2) = 0; b(n/2) = 1;
end
end