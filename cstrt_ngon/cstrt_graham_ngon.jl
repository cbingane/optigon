# 'cstrt_graham_ngon' provides the vertices coordinates (a,b) of a small
# n-gon for n = 2m and m >= 2
function cstrt_graham_ngon(n)
    if mod(n,2) == 0 && n>=4
        a = zeros(n-1); b = zeros(n-1);
        for i = 1:Int(n/2-1)
            a[i] = sin(2*pi*i/(n-1))/(2*cos(pi/(2*n-2))); a[n-i] = -a[i];
            b[i] = (1 - cos(2*pi*i/(n-1)))/(2*cos(pi/(2*n-2))); b[n-i] = b[i];
        end
        a[Int(n/2)] = 0; b[Int(n/2)] = 1;
        return(a,b)
    else
        throw(DomainError(n,"n must be an even integer >= 4"))
    end
end