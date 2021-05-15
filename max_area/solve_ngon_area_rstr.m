% 'solve_ngon_area_rstr.m' solves a convex restriction of the maximal
% area problem constructed around an n-gon (a,b)
function [xopt,yopt,Aopt] = solve_ngon_area_rstr(n,a,b)
cvx_begin
    variables x(n-1) y(n-1)
    variable u(n-2)
    maximize sum(u)
    subject to
        for i = 1:n-1
            x(i)^2 + y(i)^2 <= 1
            y(i) >= 0
        end
        % unit diameter
        for i = 2:n-1
            for j = 1:i-1
                (x(i) - x(j))^2 + (y(i) - y(j))^2 <= 1
            end
        end
        % counterclockwise order and area
        for i = 1:n-2
            u(i) >= 0
            (x(i+1) + y(i))^2 + (y(i+1) - x(i))^2 + 8*u(i) <= ...
                2*(a(i+1) - b(i))*(x(i+1) - y(i)) - (a(i+1) - b(i))^2 + ...
                2*(b(i+1) + a(i))*(y(i+1) + x(i)) - (b(i+1) + a(i))^2
        end
        % lower bound
%         sum(u) >= calc_area_ngon(a,b)
cvx_end
xopt = x; yopt = y; Aopt = calc_area_ngon(xopt,yopt);
end