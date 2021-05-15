# Small n-gon: Maximal area

model;

param n := 10;
param pi := 4*atan(1);
param s := (2*sqrt(114)-7)/22;
param t := (84*s^2-272*s+175)/(4*(22*s+7));
param c := (7792*s^4+16096*s^3+2568*s^2-6248*s+223)*pi^2/(768*(22*s+7))-(226*s^2+84*s*t-22*t^2-542*s-136*t+303)/(2*(22*s+7));
param u := s*pi/n+t*pi/n^2-c*pi/n^3;
param v := (pi/2-u)/(n/2-1);
param w := asin(sin(u)+sin(u+3*v/2)/(2*cos(v/2)))-u-v;
param Al := sin(u)+sin(2*v)-sin(v+w)+(n/2-3)*(sin(v)-tan(v/2))+(cos(v-w)-cos(2*v)-1/2)*tan(v/2);
param Au := (n*sin(pi/n)-(n-1)*tan(pi/(2*n-2)))/2;

var a {i in 1..n/2} >= 0;
var x {i in 1..n/2-1};
var y {i in 1..n/2-1};
var A {i in 1..n/2-1} >= 0;

maximize area: sum{i in 1..n/2-1} (A[i]);

subject to al: a[1] <= pi/6;
subject to aiu {i in 2..n/2}: a[i] <= pi/3;

subject to cycle: sum{i in 1..n/2} a[i] = pi/2;

subject to xi {i in 1..n/2-1}: x[i] = sum{j in 1..i} (-1)^(j-1)*sin(sum{k in 1..j} a[k]);
subject to yi {i in 1..n/2-1}: y[i] = sum{j in 1..i} (-1)^(j-1)*cos(sum{k in 1..j} a[k]);

subject to xend: x[n/2-1] = (-1)^(n/2)/2;

subject to area_l: A[1] = x[1];
subject to area_i {i in 2..n/2-2}: A[i] = x[i+1]*y[i-1]-y[i+1]*x[i-1];
subject to area_end: A[n/2-1] = -x[n/2-1]*y[n/2-2]-y[n/2-1]*x[n/2-2];

subject to bounds: Al <= sum{i in 1..n/2-1} (A[i]) <= Au;