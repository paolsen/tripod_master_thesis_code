% two_leg_semi_simple_right.
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

m_vec = [0 M 0 2*m 0 0 m 0]; % 2* because two legs on one side

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
%R = 0.12;
%R1 = R*sqrt(3);
%D1 = 0.1;
%D2 = 0.4;
%D3 = 0.4;
%M1 = 0.5;
%M2 = 1;
%G = 9.80665;
%parameter_fixated_eq = sym_expression2value(phi_eq,[g,m,M,d1,d2,r],[G,M1,M2,D1,D2,R]);
%parameter_fixated_eq = phi_eq;
theta1_s = 0;
theta2_s = pi/3;
theta3_s = pi/3;
rotationM = [cos(-pi/6), -sin(-pi/6),0;sin(-pi/6), cos(-pi/6),0;0,0,1];
%H = sqrt(dot((rotationM*LegTransform(theta1_s, theta2_s,theta3_s, D1,D2,D3,R)),([1;0;1].*(rotationM*LegTransform(theta1_s, theta2_s,theta3_s, D1,D2,D3,R)))));
%variable_fixated_eq = sym_expression2value(parameter_fixated_eq, [tha1(t),tha2(t),tha3(t),tha1_d(t),tha2_d(t),tha3_d(t),tha1_dd(t),tha2_dd(t),tha3_dd(t),thc1(t),thc3(t),thc1_d(t),thc3_d(t),thc1_dd(t),thc3_dd(t)],[theta1_s,theta2_s,theta3_s,0,0,0,0,0,0,0,0,0,0,0,0]);
% separate the ddot, dot and gravitational terms of \phi
sym H;
%H = sqrt(dot((rotationM*LegTransform(theta1_s, theta2_s,theta3_s, D1,D2,D3,R)),([1;0;1].*(rotationM*LegTransform(theta1_s, theta2_s,theta3_s, D1,D2,D3,R)))));
variable_fixated_eq = sym_expression2value(phi_eq, [tha1(t),tha2(t),tha3(t),tha1_d(t),tha2_d(t),tha3_d(t),tha1_dd(t),tha2_dd(t),tha3_dd(t),thc1(t),thc3(t),thc1_d(t),thc3_d(t),thc1_dd(t),thc3_dd(t)],[theta1_s,theta2_s,theta3_s,0,0,0,0,0,0,0,0,0,0,0,0]);


Alpha = variable_fixated_eq - sym_expression2value(variable_fixated_eq,[phi0_dd(t)],[0]);
Beta = variable_fixated_eq - sym_expression2value(variable_fixated_eq,[phi0_d(t)],[0]);
Gamma = variable_fixated_eq - Alpha - Beta;
A = simplify(Alpha/phi0_dd);
syms z
k1 = -pi/3;
k2 = pi/3;
phi1 = 0;
phi2 = pi/12;
k3 = sqrt(2)/(sqrt(pi))*(1.17679800*sin(pi/12)+0.294199*sqrt(2)*sin(pi/12) - 0.2941*sqrt(6)*sin(pi/12));
funct = k1*(z - phi2)+k2*(z-phi1) + k3*(z-phi1)*(z-phi2);
%funct = -z + acos((100*sqrt(3)*D2*M1*sin(z)-75*sqrt(2)*D2*M1*cos(z)-25*sqrt(6)*D2*M1*cos(z)-50*sqrt(6)*M1*R*cos(z)+2*sqrt(5163)*M1*sin(z)+400*D1*sin(z)*M1+200*M1*R*cos(z))/(200*D2*M2));
Gamma_fix = subs(Gamma,thc2,funct);
Gamma_fix = subs(Gamma_fix,phi0,z);
%Gamma_fix = sym_expression2value(Gamma_fix,[g,m,M,d1,d2,r],[G,M1,M2,D1,D2,R]);
%Phase portrait
