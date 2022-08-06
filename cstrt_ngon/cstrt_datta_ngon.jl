function cstrt_datta_ngon(n)
    if mod(log2(n),1)==0 && n>=16
        # initialization
        d = atan(tan(2*pi/n)*tan(pi/n))-asin(sin(2*pi/n)^2/sqrt(2+2*cos(2*pi/n)*cos(4*pi/n)));
        x = zeros(Int(3*n/8));
        for i = 1:Int(3*n/8)
            if mod(i,3) == 2
                x[i] = 2*(pi/n + (-1)^(i-1)*d);
            else
                x[i] = pi/n + (-1)^(i-1)*d;
            end
        end
        p = zeros(Int(3*n/8-1)); q = zeros(Int(3*n/8-1));
        p[1] = sin(x[1]); q[1] = cos(x[1]);
        for i = 2:Int(3*n/8-1)
            p[i] = p[i-1] - (-1)^i*sin(sum(x[1:i]));
            q[i] = q[i-1] - (-1)^i*cos(sum(x[1:i]));
        end
        # construction
        a = zeros(n-1); b = zeros(n-1);
        a[Int(n/2-1)] = p[1]; a[Int(n/2+1)] = -a[Int(n/2-1)];
        b[Int(n/2-1)] = q[1]; b[Int(n/2+1)] = b[Int(n/2-1)];
        a[n-1] = a[Int(n/2-1)]-sin(x[1]+x[2]/2); a[1] = -a[n-1];
        b[n-1] = b[Int(n/2-1)]-cos(x[1]+x[2]/2); b[1] = b[n-1];
        a[n-2] = p[2]; a[2] = -a[n-2];
        b[n-2] = q[2]; b[2] = b[n-2];
        for i = 1:Int(n/8-1)
            a[Int(mod(n/2-1+i*(3*n/2-2),n))] = p[3*i+1];
            a[Int(n-mod(n/2-1+i*(3*n/2-2),n))] = -a[Int(mod(n/2-1+i*(3*n/2-2),n))];
            b[Int(mod(n/2-1+i*(3*n/2-2),n))] = q[3*i+1];
            b[Int(n-mod(n/2-1+i*(3*n/2-2),n))] = b[Int(mod(n/2-1+i*(3*n/2-2),n))];
            a[Int(mod(n-1+i*(3*n/2-2),n))] = a[Int(mod(n/2-1+i*(3*n/2-2),n))]-(-1)^i*sin(sum(x[1:3*i+1])+x[3*i+2]/2);
            a[Int(n-mod(n-1+i*(3*n/2-2),n))] = -a[Int(mod(n-1+i*(3*n/2-2),n))];
            b[Int(mod(n-1+i*(3*n/2-2),n))] = b[Int(mod(n/2-1+i*(3*n/2-2),n))]-(-1)^i*cos(sum(x[1:3*i+1])+x[3*i+2]/2);
            b[Int(n-mod(n-1+i*(3*n/2-2),n))] = b[Int(mod(n-1+i*(3*n/2-2),n))];
            a[Int(mod(n-2+i*(3*n/2-2),n))] = p[3*i+2];
            a[Int(n-mod(n-2+i*(3*n/2-2),n))] = -a[Int(mod(n-2+i*(3*n/2-2),n))];
            b[Int(mod(n-2+i*(3*n/2-2),n))] = q[3*i+2];
            b[Int(n-mod(n-2+i*(3*n/2-2),n))] = b[Int(mod(n-2+i*(3*n/2-2),n))];
            a[Int(mod(n+i*(3*n/2-2),n))] = p[3*i];
            a[Int(n-mod(n+i*(3*n/2-2),n))] = -a[Int(mod(n+i*(3*n/2-2),n))];
            b[Int(mod(n+i*(3*n/2-2),n))] = q[3*i];
            b[Int(n-mod(n+i*(3*n/2-2),n))] = b[Int(mod(n+i*(3*n/2-2),n))];
        end
        a[Int(n/2)] = 0; b[Int(n/2)] = 1;
        return(a,b)
    else
        throw(DomainError(n,"n must be a power of 2 >= 16"))
    end
end