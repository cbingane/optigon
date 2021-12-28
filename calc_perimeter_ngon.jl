function calc_perimeter_ngon(a,b)
if length(a) == length(b) && length(a) >= 2
    L = 0;
    L = L + abs(a[1]+im*b[1]);
    for i = 2:length(a)
        L = L + abs((a[i]-a[i-1])+im*(b[i]-b[i-1]));
    end
    L = L + abs(a[length(a)]+im*b[length(a)]);
    return L
end
end