# 'cstrt_tamvakis_ngon' provides the vertices coordinates (a,b) of the
# Tamvakis n-gon Tn
# Please, cite...
#   N. K. Tamvakis. On the perimeter and the area of the convex polygon of
#   a given diameter. Bull. Greek Math. Soc, 28: 115-132, 1987.
#   C. Bingane. Tight bounds on the maximal perimeter and the maximal width
#   of convex small polygons. Journal of Global Optimization, 2022.
function cstrt_tamvakis_ngon(n)
if n>=3
    a = zeros(n-1); b = zeros(n-1);
    if mod(n,3) == 0
        a[Int(n/3)] = 1/2; a[Int(2*n/3)] = -a[Int(n/3)];
        b[Int(n/3)] = sqrt(3)/2; b[Int(2*n/3)] = b[Int(n/3)];
        if n>3
            for j = 1:Int((n-3)/3)
                a[Int(n/3+j)] = cos(pi/3+j*pi/n);
                b[Int(n/3+j)] = sin(pi/3+j*pi/n);
                a[Int(2*n/3+j)] = a[Int(n/3)] + cos(pi+j*pi/n);
                b[Int(2*n/3+j)] = b[Int(n/3)] + sin(pi+j*pi/n);
                a[Int(n/3-j)] = -a[Int(2*n/3+j)];
                b[Int(n/3-j)] = b[Int(2*n/3+j)];
            end
        end
    elseif mod(n,3) == 1
        a[Int((n-1)/3)] = 1/2; a[Int((2*n+1)/3)] = -a[Int((n-1)/3)];
        b[Int((n-1)/3)] = sqrt(3)/2; b[Int((2*n+1)/3)] = b[Int((n-1)/3)];
        for j = 1:Int((n-1)/3)
            a[Int((n-1)/3+j)] = cos(pi/3+j*pi/(n+2));
            b[Int((n-1)/3+j)] = sin(pi/3+j*pi/(n+2));
        end
        if n>4
            for j = 1:Int((n-4)/3)
                a[Int((2*n+1)/3+j)] = a[Int((n-1)/3)] + cos(pi+j*pi/(n-1));
                b[Int((2*n+1)/3+j)] = b[Int((n-1)/3)] + sin(pi+j*pi/(n-1));
                a[Int((n-1)/3-j)] = -a[Int((2*n+1)/3+j)];
                b[Int((n-1)/3-j)] = b[Int((2*n+1)/3+j)];
            end
        end
    else
        a[Int((n+1)/3)] = 1/2; a[Int((2*n-1)/3)] = -a[Int((n+1)/3)];
        b[Int((n+1)/3)] = sqrt(3)/2; b[Int((2*n-1)/3)] = b[Int((n+1)/3)];
        for j = 1:Int((n-2)/3)
            a[Int((2*n-1)/3+j)] = a[Int((n+1)/3)] + cos(pi+j*pi/(n+1));
            b[Int((2*n-1)/3+j)] = b[Int((n+1)/3)] + sin(pi+j*pi/(n+1));
            a[Int((n+1)/3-j)] = -a[Int((2*n-1)/3+j)];
            b[Int((n+1)/3-j)] = b[Int((2*n-1)/3+j)];
        end
        if n>5
            for j = 1:Int((n-5)/3)
                a[Int((n+1)/3+j)] = cos(pi/3+j*pi/(n-2));
                b[Int((n+1)/3+j)] = sin(pi/3+j*pi/(n-2));
            end
        end
    end
    return(a,b)
end
end