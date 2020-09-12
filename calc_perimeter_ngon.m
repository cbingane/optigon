function L = calc_perimeter_ngon(a,b)
if length(a) == length(b) && length(a) >= 2
    L = 0;
    L = L + abs(a(1)+1i*b(1));
    for i = 2:length(a)
        L = L + abs((a(i)-a(i-1))+1i*(b(i)-b(i-1)));
    end
    L = L + abs(a(length(a))+1i*b(length(a)));
end
end