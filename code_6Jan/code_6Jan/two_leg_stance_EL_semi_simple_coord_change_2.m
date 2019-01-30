% coord-change 2

syms phi0(t) l0(t) t
syms tha1(t) tha2(t) tha3(t)
syms thc1(t) thc2(t) thc3(t)

q = [phi0 ;l0; tha1; tha2; tha3; thc1; thc2; thc3];

syms phi0_d(t) l0_d(t)
syms tha1_d(t) tha2_d(t) tha3_d(t)
syms thc1_d(t) thc2_d(t) thc3_d(t)

q_d = [phi0_d; l0_d; tha1_d; tha2_d; tha3_d; thc1_d; thc2_d; thc3_d];

syms phi0_dd(t) l0_dd(t)
syms tha1_dd(t) tha2_dd(t) tha3_dd(t)
syms thc1_dd(t) thc2_dd(t) thc3_dd(t)

q_dd = [phi0_dd; l0_dd; tha1_dd; tha2_dd; tha3_dd; thc1_dd; thc2_dd; thc3_dd];

syms m M I
syms d1 d2 r
PI = sym('pi');

%fixated rotations in z and x
A0 = DH(0,phi0,0,0);
A1a = DH(r,0,l0,PI/12); %good
A2a = DH(0,-PI/2,d1,tha1);
A3a = DH(0,0,d2,tha2);
A4a = DH(0,0,d2,tha3);
A1c = DH(r,0,l0,-PI/2); %good
A2c = DH(0,-PI/2,d1,thc1);
A3c = DH(0,0,d2,thc2);
A4c = DH(0,0,d2,thc3);

cm_0 = four2three(A0*[0;0;0;1]);
cm_p = four2three(A0*A1a*[-r;0;0;1]);
cm_1a = four2three(A0*A1a*A2a*[-d1/2;0;0;1]);         % mass is 0 anyways
cm_2a = four2three(A0*A1a*A2a*A3a*[-d2;0;0;1]);       % knee
cm_3a = four2three(A0*A1a*A2a*A3a*A4a*[-d2/2;0;0;1]); % mass is 0 anyways
cm_1c = four2three(A0*A1c*A2c*[-d1/2;0;0;1]);         % mass is 0 anyways
cm_2c = four2three(A0*A1c*A2c*A3c*[-d2;0;0;1]);       % knee
cm_3c = four2three(A0*A1c*A2c*A3c*A4c*[-d2/2;0;0;1]); % mass is 0 anyways

m_vec = [0 M 0 m 0 0 m 0];

j1 = linVelocity(pos(A0),q);
j2 = linVelocity(pos(A0*A1a),q);
j3 = linVelocity(pos(A0*A1a*A2a),q);
j4 = linVelocity(pos(A0*A1a*A2a*A3a),q);
j5 = linVelocity(pos(A0*A1a*A2a*A3a*A4a),q);
%j6 = linVelocity(pos(A0*A1c),q); %all ready taken
j6 = linVelocity(pos(A0*A1c*A2c),q);
j7 = linVelocity(pos(A0*A1c*A2c*A3c),q);
j8 = linVelocity(pos(A0*A1c*A2c*A3c*A4c),q);

Jv = [j1 j2 j3 j4 j5 j6 j7 j8];

z_rot = sym([0;0;1]);
x_rot = sym([1;0;0]);
prism = sym([0;0;0]);
jw1 = v2zM(x_rot,9);
jw2 = v2zM([x_rot,prism],8);
jw3 = v2zM([x_rot,prism,z_rot],8);
jw4 = v2zM([x_rot,prism,z_rot,z_rot],8);
jw5 = v2zM([x_rot,prism,z_rot,z_rot,z_rot],8);
jw6 = v2zM([x_rot,prism,z_rot,z_rot,z_rot,z_rot],8);
jw7 = v2zM([x_rot,prism,z_rot,z_rot,z_rot,z_rot,z_rot],8);
jw8 = v2zM([x_rot,prism,z_rot,z_rot,z_rot,z_rot,z_rot,z_rot],8);
%jw9 = v2zM([x_rot,prism,z_rot,z_rot,z_rot,z_rot,z_rot,z_rot,z_rot],8);

%Jw = [jw1 jw2 jw3 jw4 jw5 jw6 jw7 jw8 jw9];
Jw = [jw1 jw2 jw3 jw4 jw5 jw6 jw7 jw8];

Rotations = [rot(GenRot(0,0,phi0)),rot(GenRot(sym(0),0,0)),rot(GenRot(tha1,0,0)),rot(GenRot(tha2,0,0)),rot(GenRot(tha3,0,0)),rot(GenRot(thc1,0,0)),rot(GenRot(thc2,0,0)),rot(GenRot(thc3,0,0))];
Inertia = [sym(zeros(3)), I*eye(3), sym(zeros(3)), sym(zeros(3)), sym(zeros(3)), sym(zeros(3)), sym(zeros(3)), sym(zeros(3)), sym(zeros(3))];

D = InertiaMatrix(m_vec,Jv,Jw,Rotations,Inertia);

K = 0.5*q_d.'*D*q_d;


% Potential Energy
syms g;

P = g*m_vec*[cm_0(3,1);cm_p(3,1);cm_1a(3,1);cm_2a(3,1);cm_3a(3,1);cm_1c(3,1);cm_2c(3,1);cm_3c(3,1)];

L = K - P;

eq = DeriveEL(L,q,q_d,q_dd,t);
phi_eq = eq(1);
% numerical values
R = 0.12;
%R1 = R*sqrt(3);
D1 = 0.1;
D2 = 0.4;
D3 = 0.4;
M1 = 0.1;
M2 = 10;
G = 9.80665;
parameter_fixated_eq = sym_expression2value(phi_eq,[g,m,M,d1,d2,r],[G,M1,M2,D1,D2,R]);
%parameter_fixated_eq = phi_eq;
theta1_s = 0;
theta2_s = pi/3;
theta3_s = pi/3;
rotationM = [cos(-pi/6), -sin(-pi/6),0;sin(-pi/6), cos(-pi/6),0;0,0,1];
H = sqrt(dot((rotationM*LegTransform(theta1_s, theta2_s,theta3_s, D1,D2,D3,R)),([1;0;1].*(rotationM*LegTransform(theta1_s, theta2_s,theta3_s, D1,D2,D3,R)))));
variable_fixated_eq = sym_expression2value(parameter_fixated_eq, [l0(t),l0_d(t),l0_dd(t),tha1(t),tha2(t),tha3(t),tha1_d(t),tha2_d(t),tha3_d(t),tha1_dd(t),tha2_dd(t),tha3_dd(t),thc1(t),thc3(t),thc1_d(t),thc3_d(t),thc1_dd(t),thc3_dd(t)],[H,0,0,theta1_s,theta2_s,theta3_s,0,0,0,0,0,0,0,0,0,0,0,0]);
% separate the ddot, dot and gravitational terms of \phi

Alpha = variable_fixated_eq - sym_expression2value(variable_fixated_eq,[phi0_dd(t)],[0]);
Beta = variable_fixated_eq - sym_expression2value(variable_fixated_eq,[phi0_d(t)],[0]);
Gamma = variable_fixated_eq - Alpha - Beta;
A = simplify(Alpha/phi0_dd);
syms z
k1 = -pi/3;
k2 = pi/3;
phi1 = -pi/12;
phi2 = pi/12;
k3 = sqrt(2)/(sqrt(pi))*(1.17679800*sin(pi/12)+0.294199*sqrt(2)*sin(pi/12) - 0.2941*sqrt(6)*sin(pi/12));
funct =- k1*(phi2 - z)+k2*(phi1- z) + k3*(phi1- z)*(phi2- z);
%funct = -z + acos((100*sqrt(3)*D2*M1*sin(z)-75*sqrt(2)*D2*M1*cos(z)-25*sqrt(6)*D2*M1*cos(z)-50*sqrt(6)*M1*R*cos(z)+2*sqrt(5163)*M1*sin(z)+400*D1*sin(z)*M1+200*M1*R*cos(z))/(200*D2*M2));
Gamma_fix = subs(Gamma,thc2,funct);
Gamma_fix = subs(Gamma_fix,phi0,z);
%Gamma_fix = sym_expression2value(Gamma_fix,[g,m,M,d1,d2,r],[G,M1,M2,D1,D2,R]);
%Phase portrait
f = @(t,Phi0) [Phi0(2); -1/A.*subs(Gamma_fix,z,Phi0(1))];
%f = @(t,Y) [Y(2); 1/A*sym_expression2value(Gamma_fix,[phi0(t)],[Y(1)])];
y1 = linspace(-pi/12,pi/12,20);
y2 = linspace(-0.2,0.2,20);

% creates two matrices one for all the x-values on the grid, and one for
% all the y-values on the grid. Note that x and y are matrices of the same
% size and shape, in this case 20 rows and 20 columns
[x,y] = meshgrid(y1,y2);
u = zeros(size(x));
v = zeros(size(x));

% we can use a single loop over each element to compute the derivatives at
% each point (y1, y2)
T=0; % we want the derivatives at each point at t=0, i.e. the starting time
for i = 1:numel(x)
    Yprime = f(T,[x(i); y(i)]);
    Yprime = formula(Yprime);
    u(i) = Yprime(1);
    v(i) = Yprime(2);
end
for i = 1:numel(x)
    Vmod = sqrt(u(i)^2 + v(i)^2)*0.01;
    u(i) = u(i)/Vmod;
    v(i) = v(i)/Vmod;
end
figure(1)
quiver(x,y,u,v,'color',[0,0,0]+0.7); figure(gcf)
xticks([-pi/12 -pi/24  0 pi/24  pi/12]);
xticklabels({'-\pi/12','-\pi/24','0','\pi/24','\pi/12'});
xlabel('Angle: \phi_0')
ylab = ylabel('Angular velocity: $\dot{\phi}_0$');
set(ylab, 'Interpreter', 'latex');
axis tight equal;
axis([-pi*0.01, pi*0.03, -0.1, 0.1]);

hold('on');

%forward euler solution
t = linspace(0,5,10000);
theta = zeros(1,length(t));
gamma = subs(Gamma_fix,z,theta);
phi  = zeros(1,length(t));
phi(1) = 0;
phi_d  = zeros(1,length(t));
phi_d(1) = 0.0;
phi_dd  = zeros(1,length(t));
h = (t(10000)-t(1))/10000;
t_1 = 0.5;
for i=1:(length(t)-1)
    phi_dd(i) = -1/A.*subs(Gamma_fix,z,phi(i));
    phi_d(i+1) = phi_d(i) + phi_dd(i)*h;
    phi(i+1) = phi(i) + phi_d(i)*h;
end



plot(phi,phi_d);
xlabel('\phi_0')
yl = ylabel('Angular velocity: $\dot{\phi}_0$');
set(yl, 'Interpreter', 'latex');
hl = legend('Vector field','Trajectory');
set(hl, 'Interpreter', 'latex');

figure(2)
plot(t,phi);
%hold('on');
%plot(t,0.25*pi*ones(1,length(t)),'color',[0,0,0]+0.7);
%plot(t,-0.25*pi*ones(1,length(t)),'color',[0,0,0]+0.7);
%plot(t,zeros(1,length(t)),'color',[0,0,0]);
yticks([-pi -pi/2 0 pi/2 pi]);
yticklabels({'-\pi','-\pi/2','0','\pi/2','\pi'});
axis([0 5 -pi pi]);
xlabel('time')
yl = ylabel('Angle: $\phi_0$');
set(yl, 'Interpreter', 'latex');
hl = legend('$\phi_0$');
set(hl, 'Interpreter', 'latex');


figure(3)
plot(t,subs(funct,z,phi));
hold('on');
plot(t,ones(length(t))*(-pi),'--', 'color',[0 0 0]);
xlabel('time')
y3 = ylabel('Angle: $\theta_{2c}$');
set(y3, 'Interpreter', 'latex');
%h3 = legend('$\theta_{2c}$');
%set(h3, 'Interpreter', 'latex');
