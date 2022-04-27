# 'jump_max_area_xy_sqt' solves the maximal area problem sequentially from an
# initial n-gon (a,b)
# Please, cite...
#   C. Bingane. Largest small polygons: A sequential convex optimization
#   approach. Optimization Letters, 2022.
function jump_max_area_xy_sqt(n)
if mod(n,2) == 0 && n >= 6
    (a,b) = cstrt_graham_ngon(n);
    A0 = calc_area_ngon(a,b);
    (x,y,ct) = jump_max_area_xy_rstr_quad(n,a,b);
    A = calc_area_ngon(x,y);
    tct = ct;
    k = 1;
    while (norm((x-a)+im*(y-b))/norm(x+im*y) > 1e-5)
        a = x; b = y; A0 = A;
        (x,y,ct) = jump_max_area_xy_rstr_quad(n,a,b);
        A = calc_area_ngon(x,y);
        tct = tct + ct;
        k = k+1;
    end
    ez = norm((x-a)+im*(y-b))/norm(x+im*y); eA = (A-A0)/A;
    return(A,x,y,tct,k,ez,eA)
end
end