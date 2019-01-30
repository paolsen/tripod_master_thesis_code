%two_leg_stance_semi-simple

%three-pedal model, including feet angle

syms al_0(t) th_1(t) al_1(t) al_2(t) al_3(t) al_4(t)    % initial leg (a)
syms al_5(t) al_6(t) al_7(t)                            % leg c 

q = [al_0; al_1; th_1; al_2; al_3; al_4; al_5; al_6; al_7];

syms d1 d2 d3                                       % link length
syms d                                              % platform side length
PI = sym('pi');

A0 = DH(0,al_0,0,0);                                % passive degree of freedom
A1 = DH(0,al_1,0,PI/2);                             % angle on floor
A2 = DH(0,al_2,d2,th_1 - PI/2);                     % becomes equal product of other angles
A3 = DH(0,al_3,d2,0);                               % knee to thigh                             
A4 = DH(0,al_4+PI/12,d1,PI/2);                      % hip i attatched to platform
A4_CM = DH(0,al_4,d1,PI/2);                         % hip i attatched to platform pointing towards CM


A5 = DH(0,al_5 + PI/12,d,0);                        % hip c
A6 = DH(0,al_6,d1,-PI/2);                           % hip 1 c
A7 = DH(0,al_7,d2,0);                               % hip 2 c 


% mass of links
%syms M m1 m2 m3
syms M m


% possitions of center of mass
% for simplicity centers of mass will be placed in the middle of each link.

cm_a3 = four2three(A0*A1)*[0; 0; d2;1];                           % knee
cm_a2 = four2three(A0*A1*A2)*[0; 0; d2/2;1];                      % middle of thigh
cm_a1 = four2three(A0*A1*A2*A3)*[0;0;d1/2;1];                     % hip joint mid

cm_p  = four2three(A0*A1*A2*A3*A4_CM)*[0;0;d/sqrt(3);1];          % mid platform

cm_c1 = four2three(A0*A1*A2*A3*A4*A5)*[0;0;d1/2;1];               % hip leg c
cm_c2 = four2three(A0*A1*A2*A3*A4*A5*A6)*[0;0;d2/2;1];            % mid thigh
cm_c3 = four2three(A0*A1*A2*A3*A4*A5*A6*A7) * [0;0;0;1];          % knee

% Kinetic and potential energy can be calculated
syms al_0_d(t) th_1_d(t) al_1_d(t) al_2_d(t) al_3_d(t) al_4_d(t)    % initial leg (a)
syms al_5_d(t) al_6_d(t) al_7_d(t)                                  % leg c 

q_d = [al_0_d; al_1_d; th_1_d; al_2_d; al_3_d; al_4_d; al_5_d; al_6_d; al_7_d];


j0 = linVelocity(pos(A0),[al_0; al_1; th_1; al_2; al_3; al_4; al_5; al_6; al_7]);
j1 = linVelocity(pos(A0*A1),[al_0; al_1; th_1; al_2; al_3; al_4; al_5; al_6; al_7]);
j2 = linVelocity(pos(A0*A1*A2),[al_0; al_1; th_1; al_2; al_3; al_4; al_5; al_6; al_7]);
j3 = linVelocity(pos(A0*A1*A2*A3),[al_0; al_1; th_1; al_2; al_3; al_4; al_5; al_6; al_7]);
j4 = linVelocity(pos(A0*A1*A2*A3*A4_CM),[al_0; al_1; th_1; al_2; al_3; al_4; al_5; al_6; al_7]);
j5 = linVelocity(pos(A0*A1*A2*A3*A4*A5),[al_0; al_1; th_1; al_2; al_3; al_4; al_5; al_6; al_7]);
j6 = linVelocity(pos(A0*A1*A2*A3*A4*A5*A6),[al_0; al_1; th_1; al_2; al_3; al_4; al_5; al_6; al_7]);
j7 = linVelocity(pos(A0*A1*A2*A3*A4*A5*A6*A7),[al_0; al_1; th_1; al_2; al_3; al_4; al_5; al_6; al_7]);
% linear velocity jacobian

Jv = [j0 j0 j1 j2 j3 j4 j5 j6 j7];
%Jv = [j2 j3 j4 j5 j6 j7];

% development of angular velocity Jacobian
z_rot = sym([0;0;1]);
x_rot = sym([1;0;0]);
jw1 = v2zM(x_rot,9);
jw2 = v2zM([x_rot,x_rot],9);
jw3 = v2zM([x_rot,x_rot,z_rot],9);
jw4 = v2zM([x_rot,x_rot,z_rot,x_rot],9);
jw5 = v2zM([x_rot,x_rot,z_rot,x_rot,x_rot],9);
jw6 = v2zM([x_rot,x_rot,z_rot,x_rot,x_rot,x_rot],9);
jw7 = v2zM([x_rot,x_rot,z_rot,x_rot,x_rot,x_rot,x_rot],9);
jw8 = v2zM([x_rot,x_rot,z_rot,x_rot,x_rot,x_rot,x_rot,x_rot],9);
jw9 = v2zM([x_rot,x_rot,z_rot,x_rot,x_rot,x_rot,x_rot,x_rot,x_rot],9);
Jw = [jw1,jw2,jw3,jw4,jw5,jw6,jw7,jw8,jw9];
m_vec = [0, 0, m,0,0,M,0,0,m];% modified
% Kinetic energy
Rotations = [rot(GenRot(0,0,al_0)),rot(GenRot(0,0,al_1)),rot(GenRot(th_1,0,0)),rot(GenRot(0,0,al_2)),rot(GenRot(0,0,al_3)),rot(GenRot(0,0,al_4)),rot(GenRot(0,0,al_5)),rot(GenRot(0,0,al_6)),rot(GenRot(0,0,al_7))];
%syms I i1 i2 i3
syms I
i1 = 0;
i2 = 0;
i3 = 0;
Inertia = [sym(zeros(3)),sym(zeros(3)), i3 * eye(3), i2* eye(3), i1*eye(3), I*eye(3),  i1*eye(3), i2*eye(3), i3*eye(3)];
% InertiaMatrix(massVector,LinVelJac,AngVelJac,Rotations,InertiaM)
D = InertiaMatrix(m_vec,Jv,Jw,Rotations,Inertia);

K = 0.5*q_d.'*D*q_d;

% Potential Energy
syms g;
% mass vector from first link to last
 

P = g*m_vec*[0;0;cm_a3(3,1);cm_a2(3,1);cm_a1(3,1);cm_p(3,1);cm_c1(3,1);cm_c2(3,1);cm_c3(3,1)];

L = K - P;

syms al_0_dd(t) th_1_dd(t) al_1_dd(t) al_2_dd(t) al_3_dd(t) al_4_dd(t)    % initial leg (a)
syms al_5_dd(t) al_6_dd(t) al_7_dd(t)                                  % leg c 

q_dd = [al_0_dd; al_1_dd; th_1_dd; al_2_dd; al_3_dd; al_4_dd; al_5_dd; al_6_dd; al_7_dd];

eq = DeriveEL(L,q,q_d,q_dd,t);


syms a0 a1 a2 t1 a3 a4 a5 a6 a7
syms a0d a1d a2d t1d a3d a4d a5d a6d a7d
syms a0dd a1dd a2dd t1dd a3dd a4dd a5dd a6dd a7dd
eq_nt = sym_expression2value(eq,[al_0_dd(t); al_1_dd(t); th_1_dd(t); al_2_dd(t); al_3_dd(t); al_4_dd(t); al_5_dd(t); al_6_dd(t); al_7_dd(t)] , [a0dd a1dd a2dd t1dd a3dd a4dd a5dd a6dd a7dd]);
eq_nt = sym_expression2value(eq_nt, [al_0_d(t); al_1_d(t); th_1_d(t); al_2_d(t); al_3_d(t); al_4_d(t); al_5_d(t); al_6_d(t); al_7_d(t)], [a0d a1d a2d t1d a3d a4d a5d a6d a7d]);
eq_nt = sym_expression2value(eq_nt, [al_0(t); al_1(t); th_1(t); al_2(t); al_3(t); al_4(t); al_5(t); al_6(t); al_7(t)], [a0 a1 a2 t1 a3 a4 a5 a6 a7]);
%Eq = char(eq_nt);
%fileID = fopen('equations_semi_simple.txt','w');
%fprintf(fileID,Eq);
%fclose(fileID);
D_0 = equationsToMatrix(eq_nt(1), [a0dd a1dd a2dd t1dd a3dd a4dd a5dd a6dd a7dd]);

eq_nt_1_a0 = eq_nt(1) - sym_expression2value(eq_nt(1), [a0dd; a0d; a0; cos(a0)], [ 0 0 0 0]);



