# Convex equilateral small n-gon: Maximal width

model;

param n := 8;
param pi := 4*atan(1);

var x {i in 1..n-1};
var y {i in 1..n-1};
var c;
var h {i in 1..n};
var w;

maximize width: w;

subject to y_i {i in 1..n-1}: y[i] >= 0;
subject to d_i {i in 1..n-1}: x[i]^2 + y[i]^2 <= 1;

subject to d_ij {i in 2..n-1,j in 1..i-1}: (x[i]-x[j])^2 + (y[i]-y[j])^2 <= 1;

subject to counterclockwise_order {i in 1..n-2}: x[i]*y[i+1] - y[i]*x[i+1] >= 0;

subject to side_1: c^2 = x[1]^2 + y[1]^2;
subject to side_i {i in 2..n-1}: c^2 = (x[i]-x[i-1])^2 + (y[i]-y[i-1])^2;
subject to side_n: c^2 = x[n-1]^2 + y[n-1]^2;

subject to length: 0.384462 <= c <= 0.386952;

subject to height_1: c*h[1] = -x[5]*(y[1]-y[5]) + y[5]*(x[1]-x[5]);
subject to height_2: c*h[2] = (x[1]-x[5])*(y[2]-y[5]) - (y[1]-y[5])*(x[2]-x[5]);
subject to height_3: c*h[3] = (x[2]-x[6])*(y[3]-y[6]) - (y[2]-y[6])*(x[3]-x[6]);
subject to height_4: c*h[4] = x[3]*y[4] - y[3]*x[4];
subject to height_5: c*h[5] = x[4]*y[5] - y[4]*x[5];
subject to height_6: c*h[6] = (x[5]-x[2])*(y[6]-y[2]) - (y[5]-y[2])*(x[6]-x[2]);
subject to height_7: c*h[7] = (x[6]-x[3])*(y[7]-y[3]) - (y[6]-y[3])*(x[7]-x[3]);
subject to height_8: c*h[8] = -(x[7]-x[3])*y[3] + (y[7]-y[3])*x[3];

subject to height_min {i in 1..n}: w <= h[i];

subject to bounds: w >= 0.953776;

subject to xm: x[4] = 0;

#subject to unicity: (x[1]-0.320810)^2+(y[1]-0.214)^2+(x[2]-1/2)^2+(y[2]-0.555477)^2+(x[3]-0.38411)^2+(y[3]-0.923288)^2+(y[4]-0.957567)^2+(x[5]+0.38411)^2+(y[5]-0.923288)^2+(x[6]+1/2)^2+(y[6]-0.555477)^2+(x[7]+0.320810)^2+(y[7]-0.214)^2 >= (1.5*10^(-4))^2;
