# Convex small n-gon Qn (n >= 4): Maximal perimeter

model;

param n := 8;
param pi := 4*atan(1);
param d := pi/4-asin(cos(pi/n)/sqrt(2));
param Ll := 2*n*sin(pi/(2*n))*cos(d/2);
param Lu := 2*n*sin(pi/(2*n));

var a {i in 1..n/2} >= 0;

maximize perimeter: sum{i in 1..n/2} (4*sin(a[i]/2));

subject to xu1: a[1] <= pi/6;
subject to xu {i in 1..n/2}: a[i] <= pi/3;

subject to cycle: sum{i in 1..n/2} a[i] = pi/2;

subject to x2: sum{i in 1..n/2-1} (-1)^(i-1)*sin(sum{j in 1..i} (a[j])) = (-1)^(n/2)/2;

subject to bounds: Ll <= sum{i in 1..n/2} (4*sin(a[i]/2)) <= Lu;