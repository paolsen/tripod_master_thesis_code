syms g m x x_d x_dd f my

X_DD = -m*g*10+my*x_d^4;

[q,q_d,q_dd] = ForwardEulerOrder2D1(X_DD,x,x_d,f,zeros(1,1000),[g;m;my],[9.81;1;0.001],10,10,0,5,1000);
[q_1,q_1d,q_1dd] = ForwardEulerOrder2D1(X_DD,x,x_d,f,zeros(1,1000),[g;m;my],[9.81;1;0.0001],10,-10,0,5,1000);
t = linspace(0,10,1000);
plot(t,q)
hold('on');
plot(t,q_1)