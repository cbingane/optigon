% 'solve_ngon_perimeter.m' solves sequentially the maximal perimeter
% problem from an initial n-gon (a,b)
function [x,y,L,k,ez,eL] = solve_ngon_perimeter(n)
if mod(log(n),log(2)) == 0 && n >= 4
[a,b] = cstrt_datta_ngon(n);
L0 = calc_perimeter_ngon(a,b);
xi = a; yi = b; Li = L0;
figure
hold on
plot([0;a;0],[0;b;0])
[x,y,L] = solve_ngon_perimeter_conic_restr(n,a,b);
k = 1;
plot([0;x;0],[0;y;0],'--')
while (norm((x-a)+1i*(y-b))/norm(x+1i*y) > 1e-6)
    a = x; b = y; L0 = L;
    [x,y,L] = solve_ngon_perimeter_conic_restr(n,a,b);
    k = k+1;
end
ez = norm((x-a)+1i*(y-b))/norm(x+1i*y); eL = (L-L0)/L;
sol = [(0:n-1)',[0 0; xi yi],[0 0; x y]];
plot([0;x;0],[0;y;0],'LineWidth',2)
xlabel('x')
ylabel('y')
xlim([-0.5 0.5])
ylim([0 1])
% FILE
fileID = fopen(strcat('perimeter_', num2str(n),'-gon_',datestr(now,'yyyymmddHHMMSS'),'.txt'),'w');
fprintf(fileID,'%20s\n', datetime('now'));
fprintf(fileID,'CONVEX SMALL %d-GON: MAXIMIZING PERIMETER\n',n);
fprintf(fileID,'OPTIMAL VALUES\n');
fprintf(fileID,'%13s: %15.10f\n','Initial Value',Li);
fprintf(fileID,'%13s: %15.10f\n','Optimal Value',L);
fprintf(fileID,'%13s: %15d\n','# Iterations',k);
fprintf(fileID,'OPTIMAL SOLUTIONS\n');
fprintf(fileID,'%6s %15s %15s %15s %15s\n','vertex','xi','yi','x','y');
fprintf(fileID,'%6d %15.6f %15.6f %15.6f %15.6f\n',sol');
fclose(fileID);
end
end