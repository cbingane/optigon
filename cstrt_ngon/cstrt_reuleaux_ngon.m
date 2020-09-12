function [a,b] = cstrt_reuleaux_ngon(m,k)
if mod(m,2) == 1 && m>=3 && k>=1
    a = zeros(m*k-1,1); b = a;
    if m>1
        for i = 1:m-1
            a(i*k) = sin(2*pi*i/m)/(2*cos(pi/(2*m)));
            b(i*k) = (1 - cos(2*pi*i/m))/(2*cos(pi/(2*m)));
        end
        for i = 1:(m-1)/2
            if k>1
                for j = 1:k-1
                    a((i+(m-1)/2)*k+j) = a(i*k) + cos((2*i+1+(m-3)/2+j/k)*(pi/m));
                    b((i+(m-1)/2)*k+j) = b(i*k) + sin((2*i+1+(m-3)/2+j/k)*(pi/m));
                    a(m*k-((i+(m-1)/2)*k+j)) = -a((i+(m-1)/2)*k+j);
                    b(m*k-((i+(m-1)/2)*k+j)) = b((i+(m-1)/2)*k+j);
                end
            end
        end
        if k>1
            for j = 1:k-1
                a(k*(m-1)/2+j) = cos((1+(m-3)/2+j/k)*(pi/m));
                b(k*(m-1)/2+j) = sin((1+(m-3)/2+j/k)*(pi/m));
            end
        end
    end
end
end