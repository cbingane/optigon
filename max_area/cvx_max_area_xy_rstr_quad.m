% 'cvx_max_area_xy_rstr_quad.m' solves a convex restriction of the maximal
% area problem constructed around an n-gon (a,b)
function [xopt,yopt] = cvx_max_area_xy_rstr_quad(n,a,b)
cvx_begin
%     cvx_solver_settings('MSK_IPAR_INTPNT_SOLVE_FORM', 'MSK_SOLVE_DUAL');
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
            (x(i)-y(i+1))^2+(y(i)+x(i+1))^2+8*u(i) <= ...
                2*(a(i)+b(i+1))*(x(i)+y(i+1))-(a(i)+b(i+1))^2 + ...
                2*(b(i)-a(i+1))*(y(i)-x(i+1))-(b(i)-a(i+1))^2
        end
cvx_end
xopt = x; yopt = y;
end