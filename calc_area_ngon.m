function A = calc_area_ngon(a,b)
if length(a) == length(b) && length(a) >= 2
    A = 0;
    for i = 1:length(a)-1
        A = A + (a(i)*b(i+1)-b(i)*a(i+1))/2;
    end
end
end