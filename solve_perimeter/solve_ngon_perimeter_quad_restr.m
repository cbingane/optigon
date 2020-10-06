% 'solve_ngon_perimeter_quad_restr.m' solves a quadratic convex restriction
% of the maximal perimeter problem constructed around an n-gon (a,b)
function [xopt,yopt,Lopt] = solve_ngon_perimeter_quad_restr(n,a,b)
cvx_begin
    variables x(n-1) y(n-1)
    variable v(n)
    maximize sum(v)
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
        % counterclockwise order
        for i = 1:n-2
            (x(i+1) + y(i))^2 + (y(i+1) - x(i))^2 <= ...
                2*(a(i+1) - b(i))*(x(i+1) - y(i)) - (a(i+1) - b(i))^2 + ...
                2*(b(i+1) + a(i))*(y(i+1) + x(i)) - (b(i+1) + a(i))^2
        end
        % convexity
        for i = 1:n-3
            (x(i) + sqrt(3)*y(i) + x(i+1) - sqrt(3)*y(i+1) - 2*x(i+2))^2 + ...
                (-sqrt(3)*x(i) + y(i) + sqrt(3)*x(i+1) + y(i+1) - 2*y(i+2))^2 <= ...
                2*(a(i) - sqrt(3)*b(i) + a(i+1) + sqrt(3)*b(i+1) - 2*a(i+2))*(x(i) - sqrt(3)*y(i) + x(i+1) + sqrt(3)*y(i+1) - 2*x(i+2)) - ...
                (a(i) - sqrt(3)*b(i) + a(i+1) + sqrt(3)*b(i+1) - 2*a(i+2))^2 + ...
                2*(sqrt(3)*a(i) + b(i) - sqrt(3)*a(i+1) + b(i+1) - 2*b(i+2))*(sqrt(3)*x(i) + y(i) - sqrt(3)*x(i+1) + y(i+1) - 2*y(i+2)) - ...
                (sqrt(3)*a(i) + b(i) - sqrt(3)*a(i+1) + b(i+1) - 2*b(i+2))^2
        end
        % side length
        for i = 1:n
            v(i) >= 0
        end
        v(1)^2 <= 2*a(1)*x(1) - a(1)^2 + 2*b(1)*y(1) - b(1)^2
        for i = 2:n-1
            v(i)^2 <= 2*(a(i) - a(i-1))*(x(i) - x(i-1)) - (a(i) - a(i-1))^2 + ...
                2*(b(i) - b(i-1))*(y(i) - y(i-1)) - (b(i) - b(i-1))^2
        end
        v(n)^2 <= 2*a(n-1)*x(n-1) - a(n-1)^2 + 2*b(n-1)*y(n-1) - b(n-1)^2
        % lower bound
        sum(v) >= calc_perimeter_ngon(a,b)
cvx_end
xopt = x; yopt = y; Lopt = calc_perimeter_ngon(xopt,yopt);
end