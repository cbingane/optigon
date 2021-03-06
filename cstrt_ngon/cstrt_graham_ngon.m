% 'cstrt_graham_ngon' provides the vertices coordinates (a,b) of a small
% n-gon for n = 2m and m >= 3
function [a,b] = cstrt_graham_ngon(n)
if mod(n,2) == 0 && n>=4
    a = zeros(n-1,1); b = a;
    for i = 1:n/2-1
        a(i) = sin(2*pi*i/(n-1))/(2*cos(pi/(2*n-2))); a(n-i) = -a(i);
        b(i) = (1 - cos(2*pi*i/(n-1)))/(2*cos(pi/(2*n-2))); b(n-i) = b(i);
    end
    a(n/2) = 0; b(n/2) = 1;
end
end