# Convex small n-gon Cn (n >= 16): Maximal perimeter

model;

param n := 128;
param pi := 4*atan(1);
param d := atan(tan(2*pi/n)*tan(pi/n))-asin(sin(2*pi/n)^2/sqrt(2+2*cos(2*pi/n)*cos(4*pi/n)));
param Ll := 2*n*sin(pi/(2*n))*cos(d/2);
param Lu := 2*n*sin(pi/(2*n));
param c {i in 1..3*n/8} = 1 + round((i mod 3)/3);

var a {i in 1..3*n/8} >= 0;

maximize perimeter: sum{i in 1..3*n/8} (4*c[i]*sin(a[i]/2));

subject to xu {i in 1..3*n/8}: c[i]*a[i] <= pi/3;

subject to cycle: sum{i in 1..3*n/8} (c[i]*a[i]) = pi/2;

subject to x4: sum{i in 1..3*n/8-1} (-1)^(i-1)*sin(sum{j in 1..i} (c[j]*a[j])) = 1/2;

subject to bounds: Ll <= sum{i in 1..3*n/8} (4*c[i]*sin(a[i]/2)) <= Lu;