# 'cstrt_regular_ngon' provides the vertices coordinates (a,b) of the
# regular small n-gon
function cstrt_regular_ngon(n)
    Pi = big(pi);
    if n >= 3
        a = Vector{BigFloat}(undef,n-1); b = Vector{BigFloat}(undef,n-1);
        if mod(n,2) == 1
            for i = 1:n-1
                a[i] = sin(2*Pi*i/n)/(2*cos(Pi/(2*n)));
                b[i] = (1 - cos(2*Pi*i/n))/(2*cos(Pi/(2*n)));
            end
        else
            for i = 1:n-1
                a[i] = sin(2*Pi*i/n)/2;
                b[i] = (1 - cos(2*Pi*i/n))/2;
            end
        end
        return(a,b)
    else
        throw(DomainError(n,"n must be an integer >= 3"))
    end
end