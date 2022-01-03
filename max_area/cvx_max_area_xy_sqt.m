% 'solve_ngon_area.m' solves the maximal area problem sequentially from an
% initial n-gon (a,b)
function [x,y,A,k,ez,eA] = solve_ngon_area(n)
if mod(n,2) == 0 && n >= 6
[a,b] = cstrt_bwinja_ngon(n);
A0 = calc_area_ngon(a,b);
xi = a; yi = b; Ai = A0;
figure
hold on
plot([0;a;0],[0;b;0])
[x,y,A] = solve_ngon_area_rstr(n,a,b);
k = 1;
while (norm((x-a)+1i*(y-b))/norm(x+1i*y) > 1e-5)
    a = x; b = y; A0 = A;
    [x,y,A] = solve_ngon_area_rstr(n,a,b);
    k = k+1;
end
ez = norm((x-a)+1i*(y-b))/norm(x+1i*y); eA = (A-A0)/A;
sol = [(0:n-1)',[0 0; xi yi],[0 0; x y]];
plot([0;x;0],[0;y;0],'LineWidth',2)
xlabel('x')
ylabel('y')
xlim([-0.5 0.5])
ylim([0 1])
% FILE
fileID = fopen(strcat('area_', num2str(n),'-gon_',datestr(now,'yyyymmddHHMMSS'),'.txt'),'w');
fprintf(fileID,'%20s\n', datetime('now'));
fprintf(fileID,'SMALL %d-GON: MAXIMIZING AREA\n',n);
fprintf(fileID,'OPTIMAL VALUES\n');
fprintf(fileID,'%13s: %15.10f\n','Initial Value',Ai);
fprintf(fileID,'%13s: %15.10f\n','Optimal Value',A);
fprintf(fileID,'%13s: %15d\n','# Iterations',k);
fprintf(fileID,'OPTIMAL SOLUTIONS\n');
fprintf(fileID,'%6s %15s %15s %15s %15s\n','vertex','xi','yi','x','y');
fprintf(fileID,'%6d %15.6f %15.6f %15.6f %15.6f\n',sol');
fclose(fileID);
end
end
