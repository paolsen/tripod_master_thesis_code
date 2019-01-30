% calculation of torques for legs
syms th1 th2 th3 t d1 d2 d3 d d0 m g h f_x f_y f_z


J = jacobian([(d1+d2*cos(th2)+d3*cos(th2+th3))*cos(th1),(d1+d2*cos(th2)+d3*cos(th2+th3))*sin(th1),-(d2*sin(th2)+d3*sin(th2+th3))],[th1,th2,th3]);
J_funct = @(th1, th2,th3) J;
n = [0; 0; m*g/2*(d/h)];
n_c = [0;0;m*g*(1-d/h)];

f_p = [-f_x;0;-f_z];

%tau_a = (n+rot(GenRot(2/3*sym(pi),0,0))*f_p).'*J;
%tau_b = (n+rot(GenRot(4/3*sym(pi),0,0))*f_p).'*J;
tau_c = (n_c+ f_p).'*J;

%tau_a.'
%tau_b.'
tau_c.';
Tc = tau_c.';

%test = sym_expression2value(tau_c,[d3, g, m, th2, th3, d(t), h],[1, 1, 1, 1, 1, 1, 1])
%test = sym_expression2value(tau_c.',[d3, g, m, th1, h],[1, 9.81, 1, 0, 1])
D3 = 0.4;
D2 = 0.4;
tau_fix_param = sym_expression2value(Tc,[m, g, h, d1, d2, d3, th1],[1,9.80665,1, 0.1,0.4,0.4,0]);
t_0 = 0;
t_1 = 2;
steps = 10000;
t = linspace(t_0,t_1,steps);
dt = (t_1 -t_0)/steps;
theta_2 = zeros(1,steps);
theta_2_0 = pi/3;
theta_3 = zeros(1,steps);
theta_3_0 = pi/3;
theta_3_1 = pi/6;
d_h = zeros(1,steps);
surface_angle = theta_2_0 + theta_3_0;
H0 = sin(theta_3_0+theta_2_0)*D3 + sin(theta_2_0)*D2;
t3 = @(t) theta_3_0/t_1.*(t_1-t)+theta_3_1/t_1.*t;
tau2 = zeros(1,steps);
tau3 = zeros(1,steps);
torque = zeros(3,steps);
for i=1:(steps-1)
    theta_3(i) = t3(t(i));
    %theta_2(i) = - theta_3(i) -atan2(-(0.064*cos(theta_3(i)) + 0.064 - sqrt(0.004096*cos(theta_3(i))^2 +0.0512*cos(theta_3(i))*sin(theta_3(i))^2 + 0.0512*sin(theta_3(i))^2 -0.004096)/(0.32*cos(theta_3(i)) + 0.32)),-1/(sin(theta_3(i)))*(2.5*(-0.4*cos(theta_3(i))*(0.064*cos(theta_3(i))+0.064 - sqrt(0.004096*cos(theta_3(i))^2 +0.0512*cos(theta_3(i))*sin(theta_3(i))^2 + 0.0512*sin(theta_3(i))^2 -0.004096)))/(0.32*cos(theta_3(i)) + 0.32) - 0.4*(0.064*cos(theta_3(i))+ +0.064 - sqrt(0.004096*cos(theta_3(i))^2 +0.0512*cos(theta_3(i))*sin(theta_3(i))^2 + 0.0512*sin(theta_3(i))^2 -0.004096)))/(0.32*cos(theta_3(i)) + 0.32) +0.16);
    %theta_2(i) = atan2((H0*sin(theta_3(i))/2 - sqrt(sin(theta_3(i))^2*H0^2 - 2*cos(theta_3(i))*sin(theta_3(i))^2 + 2*cos(theta_3(i))*H0^2+2*sin(theta_3(i))^2 - 2*H0^2)/2)/(sin(theta_3(i))), -(((H0*sin(theta_3(i))/2 - sqrt(sin(theta_3(i))^2*H0^2 - 2*cos(theta_3(i))*sin(theta_3(i))^2 + 2*cos(theta_3(i))*H0^2+2*sin(theta_3(i))^2 - 2*H0^2)/2)*cos(theta_3(i))) - (H0*sin(theta_3(i))/2 - sqrt(sin(theta_3(i))^2*H0^2 - 2*cos(theta_3(i))*sin(theta_3(i))^2 + 2*cos(theta_3(i))*H0^2+2*sin(theta_3(i))^2 - 2*H0^2)/2)))/(sin(theta_3(i))^2);
    %theta_2(i) = (pi*sin(pi/6) +6*cos(pi/6)-2*sqrt(18*sin(pi/6)^2 - 18*sin(pi/6)*H0 + 3*pi*sin(pi/6) + 18*sin(pi/6)*theta_3(i) + 9*cos(pi/6)^2 + 18*cos(pi/6) + 9)+6)/(6*sin(pi/6));
    %theta_2(i) = solve(-th2 + asin(H0/D2 - sin(th2))==theta_3(i), th2);
    %theta_2(i) = 0.1e1 / cos(pi / 0.6e1) * (-0.8e1 * sin(pi / 0.6e1) ^ 3 - 0.24e2 * cos(pi / 0.6e1) ^ 2 * H0 + 0.4e1 * cos(pi / 0.6e1) ^ 2 * pi + 0.24e2 * cos(pi / 0.6e1) ^ 2 * theta_3(i) + 0.4e1 * sqrt(-0.24e2 * sin(pi / 0.6e1) ^ 4 + 0.24e2 * sin(pi / 0.6e1) ^ 3 * H0 - 0.4e1 * sin(pi / 0.6e1) ^ 3 * pi - 0.24e2 * sin(pi / 0.6e1) ^ 3 * theta_3(i) - 0.48e2 * sin(pi / 0.6e1) ^ 2 * cos(pi / 0.6e1) ^ 2 - 0.32e2 * cos(pi / 0.6e1) ^ 4 + 0.36e2 * cos(pi / 0.6e1) ^ 2 * H0 ^ 2 - 0.12e2 * cos(pi / 0.6e1) ^ 2 * H0 * pi - 0.72e2 * cos(pi / 0.6e1) ^ 2 * H0 * theta_3(i) + cos(pi / 0.6e1) ^ 2 * pi ^ 2 + 0.12e2 * cos(pi / 0.6e1) ^ 2 * pi * theta_3(i) + 0.36e2 * cos(pi / 0.6e1) ^ 2 * theta_3(i) ^ 2 - 0.96e2 * sin(pi / 0.6e1) ^ 2 * cos(pi / 0.6e1) + 0.72e2 * sin(pi / 0.6e1) * cos(pi / 0.6e1) * H0 - 0.12e2 * cos(pi / 0.6e1) * pi * sin(pi / 0.6e1) - 0.72e2 * sin(pi / 0.6e1) * cos(pi / 0.6e1) * theta_3(i) - 0.96e2 * cos(pi / 0.6e1) ^ 3 - 0.12e2 * sin(pi / 0.6e1) ^ 2 - 0.96e2 * cos(pi / 0.6e1) ^ 2 - 0.32e2 * cos(pi / 0.6e1)) * cos(pi / 0.6e1) - 0.24e2 * cos(pi / 0.6e1) * sin(pi / 0.6e1)) ^ (0.1e1 / 0.3e1) / 0.2e1 + 0.2e1 * (sin(pi / 0.6e1) ^ 2 + 0.2e1 * cos(pi / 0.6e1) ^ 2 + 0.2e1 * cos(pi / 0.6e1)) / cos(pi / 0.6e1) * (-0.8e1 * sin(pi / 0.6e1) ^ 3 - 0.24e2 * cos(pi / 0.6e1) ^ 2 * H0 + 0.4e1 * cos(pi / 0.6e1) ^ 2 * pi + 0.24e2 * cos(pi / 0.6e1) ^ 2 * theta_3(i) + 0.4e1 * sqrt(-0.24e2 * sin(pi / 0.6e1) ^ 4 + 0.24e2 * sin(pi / 0.6e1) ^ 3 * H0 - 0.4e1 * sin(pi / 0.6e1) ^ 3 * pi - 0.24e2 * sin(pi / 0.6e1) ^ 3 * theta_3(i) - 0.48e2 * sin(pi / 0.6e1) ^ 2 * cos(pi / 0.6e1) ^ 2 - 0.32e2 * cos(pi / 0.6e1) ^ 4 + 0.36e2 * cos(pi / 0.6e1) ^ 2 * H0 ^ 2 - 0.12e2 * cos(pi / 0.6e1) ^ 2 * H0 * pi - 0.72e2 * cos(pi / 0.6e1) ^ 2 * H0 * theta_3(i) + cos(pi / 0.6e1) ^ 2 * pi ^ 2 + 0.12e2 * cos(pi / 0.6e1) ^ 2 * pi * theta_3(i) + 0.36e2 * cos(pi / 0.6e1) ^ 2 * theta_3(i) ^ 2 - 0.96e2 * sin(pi / 0.6e1) ^ 2 * cos(pi / 0.6e1) + 0.72e2 * sin(pi / 0.6e1) * cos(pi / 0.6e1) * H0 - 0.12e2 * cos(pi / 0.6e1) * pi * sin(pi / 0.6e1) - 0.72e2 * sin(pi / 0.6e1) * cos(pi / 0.6e1) * theta_3(i) - 0.96e2 * cos(pi / 0.6e1) ^ 3 - 0.12e2 * sin(pi / 0.6e1) ^ 2 - 0.96e2 * cos(pi / 0.6e1) ^ 2 - 0.32e2 * cos(pi / 0.6e1)) * cos(pi / 0.6e1) - 0.24e2 * cos(pi / 0.6e1) * sin(pi / 0.6e1)) ^ (-0.1e1 / 0.3e1) - (-cos(pi / 0.6e1) * pi + 0.6e1 * sin(pi / 0.6e1)) / cos(pi / 0.6e1) / 0.6e1;
    theta_2(i+1) = real(asin(H0/D2 - sin(theta_3(i) + theta_2(i))));
    d_h(i) = abs(0.1 + 0.12 + (D2*cos(theta_2(i))+D3*cos(theta_2(i)+theta_3(i))));
    torque = sym_expression2value(((n_c+ [0;0;0]).'*J).', [m, g, h, d1, d2, d3, th1, th2, th3,d],[0.5,9.80665,0.5, 0.1,0.4,0.4,0,theta_2(i),theta_3(i),d_h(i)]);
    tau2(i) = torque(2);
    tau3(i) = torque(3);
end
figure(1)
plot(t,-(D2*sin(theta_2)+D3*sin(theta_2+theta_3)))
%plot(t,0.1 +(D2*cos(theta_2)+D3*cos(theta_2+theta_3)))
%plot(0.1+ (D2*cos(theta_2)+D3*cos(theta_2+theta_3)),H0-(D2*sin(theta_2)+D3*sin(theta_2+theta_3)));
figure(2)
plot(t,theta_2, t,theta_3);
yticks([0 pi/6 pi/3 pi/2]);
yticklabels({'0','\pi/6','\pi/3','\pi/2'});
axis([0 2 0 pi/2]);
xlabel('time')
yl = ylabel('Angle:');
set(yl, 'Interpreter', 'latex');
hl = legend('$\theta_2$', '$\theta_3$');
set(hl, 'Interpreter', 'latex');

figure(3)
plot(t,0.1+ (D2*cos(theta_2)+D3*cos(theta_2+theta_3)));
%yticks([0 pi/6 pi/3 pi/2]);
%yticklabels({'0','\pi/6','\pi/3','\pi/2'});
axis([0 2 0 pi/2]);
xlabel('time')
yl = ylabel('Angle:');
set(yl, 'Interpreter', 'latex');
hl = legend('$d(t)$');
set(hl, 'Interpreter', 'latex');

figure(4)
plot(t,tau2);
hold('on');
plot(t,tau3);
yl = ylabel('Torque:');
set(yl, 'Interpreter', 'latex');
hl = legend('$\tau_2$', '$\tau_3$');
set(hl, 'Interpreter', 'latex');
xlabel('time')
hold('off');
