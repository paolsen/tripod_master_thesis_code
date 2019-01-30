%three-pedal model, including feet angle

syms th_a0(t) al_a0(t) th_a1(t) al_a1(t) al_a2(t) al_a3(t)   % initial leg (a)
syms al_b4(t) al_b5(t) al_b6(t)                     % leg b
syms al_c4(t) al_c5(t) al_c6(t)                     % leg c 

q = [th_a0; al_a0; th_a1; al_a1; al_a2; al_a3; al_b4; al_b5; al_b6; al_c4; al_c5; al_c6];

syms d1 d2 d3                                       % link length
syms d                                              % platform side length
PI = sym('pi');

A0_0 = DH(0,0,0,th_a0);                             % used in calculatin jacobian
A0 = DH(0,al_a0,0,th_a0);                           % angle on floor
A1_0 = DH(0,0,0,th_a1);                             % used in calculatin jacobian
A1 = DH(0,al_a1,d3,th_a1);                          % knee
A2 = DH(0,al_a2,d2,0);                              % hip 2
A3 = DH(0,al_a3,d1,PI/2);                           % hip i attatched to platform

A3b = DH(0,al_a3 - PI/6,d1,PI/2);                   % towards leg b

A3c = DH(0,al_a3 + PI/6,d1,PI/2);                   % towards leg c

A4b = DH(0,al_b4 + PI/6,d,0);                       % hip 1 b
A5b = DH(0,al_b5,d1,-PI/2);                         % hip 2 b
A6b = DH(0,al_b6,d2,0);                             % knee b

A4c = DH(0,al_c4 - PI/6,d,0);                       % hip 1 c
A5c = DH(0,al_c5,d1,-PI/2);                         % hip 2 c
A6c = DH(0,al_c6,d2,0);                             % knee c


% mass of links
syms M m1 m2 m3
% mass vector from first link to last
m_vec = [0, 0, m3,m2,m1,M,m1,m2,m3,m1,m2,m3];

% possitions of center of mass
% for simplicity centers of mass will be placed in the middle of each link.

cm_a3 = four2three(A0)*[0; 0; d3/2;1];                              % middle of leg
cm_a2 = four2three(A0*A1)*[0; 0; d2/2;1];                           % middle of thigh
cm_a1 = four2three(A0*A1*A2)*[0;0;d1/2;1];                          % hip

cm_p  = four2three(A0*A1*A2*A3)*[0;0;d/sqrt(3);1];                  % mid platform

cm_b1 = four2three(A0*A1*A2*A3b*A4b)*[0;0;d1/2;1];                  % hip leg b
cm_b2 = four2three(A0*A1*A2*A3b*A4b*A5b)*[0;0;d2/2;1];              % middel thigh
cm_b3 = four2three(A0*A1*A2*A3b*A4b*A5b*A6b) * [0;0;d3/2;1];        % middel of leg

cm_c1 = four2three(A0*A1*A2*A3c*A4c)*[0;0;d1/2;1];                  % hip leg c
cm_c2 = four2three(A0*A1*A2*A3c*A4c*A5c)*[0;0;d2/2;1];              % middel thigh
cm_c3 = four2three(A0*A1*A2*A3c*A4c*A5c*A6c) * [0;0;d3/2;1];        % middel of leg

% Kinetic and potential energy can be calculated

syms th_a0_d(t) al_a0_d(t) th_a1_d(t) al_a1_d(t) al_a2_d(t) al_a3_d(t)         % initial leg (a)
syms al_b4_d(t) al_b5_d(t) al_b6_d(t)                               % leg b
syms al_c4_d(t) al_c5_d(t) al_c6_d(t)                               % leg c

q_d = [th_a0_d; al_a0_d; al_a1_d;th_a1_d; al_a2_d; al_a3_d; al_b4_d; al_b5_d; al_b6_d; al_c4_d; al_c5_d; al_c6_d];


j0_0 = linVelocity(pos(A0_0),[th_a0;al_a0;th_a1;al_a1;al_a2;al_a3;al_b4;al_b5;al_b6;al_c4;al_c5;al_c6]);
j0 = linVelocity(pos(A0),[th_a0;al_a0;th_a1;al_a1;al_a2;al_a3;al_b4;al_b5;al_b6;al_c4;al_c5;al_c6]);
j1_0 = linVelocity(pos(A1_0),[th_a0;al_a0;th_a1;al_a1;al_a2;al_a3;al_b4;al_b5;al_b6;al_c4;al_c5;al_c6]);
j1 = linVelocity(pos(A0*A1),[th_a0;al_a0;th_a1;al_a1;al_a2;al_a3;al_b4;al_b5;al_b6;al_c4;al_c5;al_c6]);
j2 = linVelocity(pos(A0*A1*A2),[th_a0;al_a0;th_a1;al_a1;al_a2;al_a3;al_b4;al_b5;al_b6;al_c4;al_c5;al_c6]);
j3 = linVelocity(pos(A0*A1*A2*A3),[th_a0;al_a0;th_a1;al_a1;al_a2;al_a3;al_b4;al_b5;al_b6;al_c4;al_c5;al_c6]);
jb4 = linVelocity(pos(A0*A1*A2*A3b*A4b),[th_a0;al_a0;th_a1;al_a1;al_a2;al_a3;al_b4;al_b5;al_b6;al_c4;al_c5;al_c6]);
jb5 = linVelocity(pos(A0*A1*A2*A3b*A4b*A5b),[th_a0;al_a0;th_a1;al_a1;al_a2;al_a3;al_b4;al_b5;al_b6;al_c4;al_c5;al_c6]);
jb6 = linVelocity(pos(A0*A1*A2*A3b*A4b*A5b*A6b),[th_a0;al_a0;th_a1;al_a1;al_a2;al_a3;al_b4;al_b5;al_b6;al_c4;al_c5;al_c6]);
jc4 = linVelocity(pos(A0*A1*A2*A3c*A4c),[th_a0;al_a0;th_a1;al_a1;al_a2;al_a3;al_b4;al_b5;al_b6;al_c4;al_c5;al_c6]);
jc5 = linVelocity(pos(A0*A1*A2*A3c*A4c*A5c),[th_a0;al_a0;th_a1;al_a1;al_a2;al_a3;al_b4;al_b5;al_b6;al_c4;al_c5;al_c6]);
jc6 = linVelocity(pos(A0*A1*A2*A3c*A4c*A5c*A6c),[th_a0;al_a0;th_a1;al_a1;al_a2;al_a3;al_b4;al_b5;al_b6;al_c4;al_c5;al_c6]);
% linear velocity jacobian

Jv = [j0_0 j0 j1_0 j1 j2 j3 jb4 jb5 jb6 jc4 jc5 jc6]; 

% development of angular velocity Jacobian
z_rot = sym([0;0;1]);
x_rot = sym([1;0;0]);
jw1 = v2zM(z_rot,12);
jw2 = v2zM([z_rot,x_rot],12);
jw3 = v2zM([z_rot,x_rot,z_rot],12);
jw4 = v2zM([z_rot,x_rot,z_rot,x_rot],12);
jw5 = v2zM([z_rot,x_rot,z_rot,x_rot,x_rot],12);
jw6 = v2zM([z_rot,x_rot,z_rot,x_rot,x_rot,x_rot],12);
jw7 = v2zM([z_rot,x_rot,z_rot,x_rot,x_rot,x_rot,x_rot],12);
jw8 = v2zM([z_rot,x_rot,z_rot,x_rot,x_rot,x_rot,x_rot,x_rot],12);
jw9 = v2zM([z_rot,x_rot,z_rot,x_rot,x_rot,x_rot,x_rot,x_rot,x_rot],12);
jw10 = v2zM([z_rot,x_rot,z_rot,x_rot,x_rot,x_rot,x_rot,x_rot,x_rot,x_rot],12);
jw11 = v2zM([z_rot,x_rot,z_rot,x_rot,x_rot,x_rot,x_rot,x_rot,x_rot,x_rot,x_rot],12);
jw12 = v2zM([z_rot,x_rot,z_rot,x_rot,x_rot,x_rot,x_rot,x_rot,x_rot,x_rot,x_rot,x_rot],12);
Jw = [jw1,jw2,jw3,jw4,jw5,jw6,jw7,jw8,jw9,jw10,jw11,jw12];

% Kinetic energy
Rotations = [rot(GenRot(th_a0,0,0)),rot(GenRot(0,0,al_a0)),rot(GenRot(th_a1,0,0)),rot(GenRot(0,0,al_a1)),rot(GenRot(0,0,al_a2)),rot(GenRot(0,0,al_a3)),rot(GenRot(0,0,al_b4)),rot(GenRot(0,0,al_b5)),rot(GenRot(0,0,al_b6)),rot(GenRot(0,0,al_c4)),rot(GenRot(0,0,al_c5)),rot(GenRot(0,0,al_c6))];
syms I i1 i2 i3
Inertia = [sym(zeros(3)),sym(zeros(3)), i3 * eye(3), i2* eye(3), i1*eye(3), I*eye(3), i1*eye(3), i2*eye(3), i3*eye(3),  i1*eye(3), i2*eye(3), i3*eye(3)];
% InertiaMatrix(massVector,LinVelJac,AngVelJac,Rotations,InertiaM)
D = InertiaMatrix(m_vec,Jv,Jw,Rotations,Inertia);

K = 0.5*q_d.'*D*q_d;

% Potential Energy
syms g;

P = g*m_vec*[0;0;cm_a3(3,1);cm_a2(3,1);cm_a1(3,1);cm_p(3,1);cm_b1(3,1);cm_b2(3,1);cm_b3(3,1);cm_c1(3,1);cm_c2(3,1);cm_c3(3,1)];

L = K - P;

syms th_a0_dd(t) al_a0_dd(t) th_a1_dd(t) al_a1_dd(t) al_a2_dd(t) al_a3_dd(t)         % initial leg (a)
syms al_b4_dd(t) al_b5_dd(t) al_b6_dd(t)                               % leg b
syms al_c4_dd(t) al_c5_dd(t) al_c6_dd(t)

q_dd = [th_a0_dd; al_a0_dd; al_a1_dd;th_a1_dd; al_a2_dd; al_a3_dd; al_b4_dd; al_b5_dd; al_b6_dd; al_c4_dd; al_c5_dd; al_c6_dd];

eq = DeriveEL(L,q,q_d,q_dd,t);

Eq = char(eq);
fileID = fopen('equations.txt','w');
fprintf(fileID,Eq);
fclose(fileID);

