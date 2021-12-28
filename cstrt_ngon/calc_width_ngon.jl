function calc_width_ngon(a,b)
if length(a) == length(b) && length(a) >= 2
    # SIDES
    v = zeros(length(a)+1);
    v[1] = abs(a[1]+im*b[1]);
    for i = 2:length(a)
        v[i] = abs((a[i]-a[i-1])+im*(b[i]-b[i-1]));
    end
    v[length(a)+1] = abs(a[length(a)]+im*b[length(a)]);
    # HEIGHTS
    h = zeros(length(a)+1);
    # h1
    j = 1;
    while j <= length(a)
        h1 = abs(a[j]*(b[1]-b[j])-b[j]*(a[1]-a[j]))/v[1];
        if h1 > h[1]
            h[1] = h1;
        end
        j = j+1;
    end
    # hi
    for i = 2:length(a)
        h[i] = abs(a[i-1]*b[i]-b[i-1]*a[i])/v[i];
        j = 1;
        while j <= length(a)
            hi = abs((a[i-1]-a[j])*(b[i]-b[j])-(b[i-1]-b[j])*(a[i]-a[j]))/v[i];
            if hi > h[i]
                h[i] = hi;
            end
            j = j+1;
        end
    end
    # hn
    j = 1;
    while j <= length(a)
        hn = abs((a[length(a)]-a[j])*b[j]-(b[length(a)]-b[j])*a[j])/v[length(a)+1];
        if hn > h[length(a)+1]
            h[length(a)+1] = hn;
        end
        j = j+1;
    end
    # WIDTH
    W = minimum(h);
    return W
end
end