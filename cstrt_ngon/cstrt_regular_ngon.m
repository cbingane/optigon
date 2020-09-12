function [a,b] = cstrt_regular_ngon(n)
if n>=3
    a = zeros(n-1,1); b = a;
    if mod(n,2) == 1
        for i = 1:n-1
            a(i) = sin(2*pi*i/n)/(2*cos(pi/(2*n)));
            b(i) = (1 - cos(2*pi*i/n))/(2*cos(pi/(2*n)));
        end
    else
        for i = 1:n-1
            a(i) = sin(2*pi*i/n)/2;
            b(i) = (1 - cos(2*pi*i/n))/2;
        end
    end
end
end