% 'cstrt_tamvakis_ngon' provides the vertices coordinates (a,b) of the
% Tamvakis n-gon
% Please, cite...
%   N. K. Tamvakis. On the perimeter and the area of the convex polygon of
%   a given diameter. Bull. Greek Math. Soc, 28: 115-132, 1987.
function [a,b] = cstrt_tamvakis_ngon(n)
if n>=3
    a = zeros(n-1,1); b = a;
    if mod(n,3) == 0
        a(n/3) = 1/2; a(2*n/3) = -a(n/3);
        b(n/3) = sqrt(3)/2; b(2*n/3) = b(n/3);
        if n>3
            for j = 1:(n-3)/3
                a(n/3+j) = cos(pi/3+j*pi/n);
                b(n/3+j) = sin(pi/3+j*pi/n);
                a(2*n/3+j) = a(n/3) + cos(pi+j*pi/n);
                b(2*n/3+j) = b(n/3) + sin(pi+j*pi/n);
                a(n/3-j) = -a(2*n/3+j);
                b(n/3-j) = b(2*n/3+j);
            end
        end
    elseif mod(n,3) == 1
        a((n-1)/3) = 1/2; a((2*n+1)/3) = -a((n-1)/3);
        b((n-1)/3) = sqrt(3)/2; b((2*n+1)/3) = b((n-1)/3);
        for j = 1:(n-1)/3
            a((n-1)/3+j) = cos(pi/3+j*pi/(n+2));
            b((n-1)/3+j) = sin(pi/3+j*pi/(n+2));
        end
        if n>4
            for j = 1:(n-4)/3
                a((2*n+1)/3+j) = a((n-1)/3) + cos(pi+j*pi/(n-1));
                b((2*n+1)/3+j) = b((n-1)/3) + sin(pi+j*pi/(n-1));
                a((n-1)/3-j) = -a((2*n+1)/3+j);
                b((n-1)/3-j) = b((2*n+1)/3+j);
            end
        end
    else
        a((n+1)/3) = 1/2; a((2*n-1)/3) = -a((n+1)/3);
        b((n+1)/3) = sqrt(3)/2; b((2*n-1)/3) = b((n+1)/3);
        for j = 1:(n-2)/3
            a((2*n-1)/3+j) = a((n+1)/3) + cos(pi+j*pi/(n+1));
            b((2*n-1)/3+j) = b((n+1)/3) + sin(pi+j*pi/(n+1));
            a((n+1)/3-j) = -a((2*n-1)/3+j);
            b((n+1)/3-j) = b((2*n-1)/3+j);
        end
        if n>5
            for j = 1:(n-5)/3
                a((n+1)/3+j) = cos(pi/3+j*pi/(n-2));
                b((n+1)/3+j) = sin(pi/3+j*pi/(n-2));
            end
        end
    end
end
end