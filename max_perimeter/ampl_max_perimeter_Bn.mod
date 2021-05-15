# Convex small n-gon Bn (n >= 8): Maximal perimeter

model;

param n := 8;
param pi := 4*atan(1);
param d := pi/n - asin(sin(2*pi/n)/2);
param Ll := 2*n*sin(pi/(2*n))*cos(d/2);
param Lu := 2*n*sin(pi/(2*n));

var a {i in 1..n/4+1} >= 0;

maximize perimeter: 4*sin(a[1]/2) + sum{i in 2..n/4} (8*sin(a[i]/2)) + 4*sin(a[n/4+1]/2);

subject to xui {i in 1..n/4}: a[i] <= pi/6;
subject to xu3: a[n/4+1] <= pi/3;

subject to cycle: a[1] + sum{i in 2..n/4} (2*a[i]) + a[n/4+1] = pi/2;

subject to x6: sin(a[1]) + sum{i in 2..n/4} (-1)^(i-1)*sin(a[1] + sum{j in 2..i} (2*a[j])) = -1/2;

subject to bound: Ll <= 4*sin(a[1]/2) + sum{i in 2..n/4} (8*sin(a[i]/2)) + 4*sin(a[n/4+1]/2) <= Lu;