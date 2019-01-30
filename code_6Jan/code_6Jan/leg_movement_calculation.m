% leg movement calculation
syms th1(t) th2(t) th3(t) d1 d2 d3 i1 i2 i3 m1 m2 m3 g t
PI = sym('pi');

q = [th1; th2; th3];

A1 = DH(0,-PI/2,d1,th1);
A2 = DH(0,0,d2,th2);
A3 = DH(0,0,d3,th3);



cm_1 = four2three(A1*[-0.5*d1;0;0;1]);       % must check
cm_2 = four2three(A1*A2*[-0.5*d2;0;0;1]);    % must check 
cm_3 = four2three(A1*A2*A3*[-0.5*d3;0;0;1]); % must check 

syms th1_d(t) th2_d(t) th3_d(t);
q_d = [th1_d;th2_d;th3_d];

j1 = linVelocity(pos(A1),q);
j2 = linVelocity(pos(A1*A2),q);
j3 = linVelocity(pos(A1*A2*A3),q);

Jv = [j1 j2 j3];

z_rot = sym([0;0;1]);
x_rot = sym([1;0;0]);
jw1 = v2zM(z_rot,3);
jw2 = v2zM([z_rot z_rot],3);
jw3 = v2zM([z_rot z_rot z_rot],3);
Jw = [jw1,jw2,jw3];

Rotations = [rot(GenRot(th1,0,0)),rot(GenRot(th2,0,0)),rot(GenRot(th3,0,0))];

Inertia = [i1 * eye(3), i2* eye(3), i3*eye(3)];

m_vec = [m1, m2, m3];

D = InertiaMatrix(m_vec,Jv,Jw,Rotations,Inertia);

K = 0.5*q_d.'*D*q_d;


P = g*m_vec*[cm_1(3,1);cm_2(3,1);cm_3(3,1)];

L = K-P;

syms th1_dd(t) th2_dd(t) th3_dd(t);
q_dd = [th1_dd; th2_dd; th3_dd];

eq = DeriveEL(L,q,q_d,q_dd,t);
%eq = sym_expression2value(eq,[d1 d2 d3 i1 i2 i3 m1 m2 m3 g], [1 1 1 1 1 1 1 1 1 10]);
Eq = char(eq);
fileID = fopen('leg_equations_simplified.txt','w');
fprintf(fileID,Eq);
fclose(fileID);
% making the functions time independant
syms theta_1 theta_2 theta_3 theta_1_d theta_2_d theta_3_d theta_1_dd theta_2_dd theta_3_dd
eq = sym_expression2value(eq,[th1(t) th2(t) th3(t) th1_d(t) th2_d(t) th3_d(t) th1_dd(t) th2_dd(t) th3_dd(t)], [theta_1 theta_2 theta_3 theta_1_d theta_2_d theta_3_d theta_1_dd theta_2_dd theta_3_dd]);
% separating the equations in two matrices and one vector

D = equationsToMatrix(eq, [theta_1_dd theta_2_dd theta_3_dd]);
eq = sym_expression2value(eq,[theta_1_dd theta_2_dd theta_3_dd], [0 0 0]);
syms theta_1_ds theta_2_ds theta_3_ds
%eq = sym_expression2value(eq,[theta_1_d^2 theta_2_d^2 theta_3_d^2], [theta_1_ds theta_2_ds theta_3_ds]);
% here comes faulty
neg = sym_expression2value(eq,[theta_1_d theta_2_d theta_3_d],[0 0 0]);
c = eq - neg;
%C = equationsToMatrix(eq,[theta_1_ds theta_2_ds theta_3_ds]);
G = neg;
%G = sym_expression2value(eq, [theta_1_ds theta_2_ds theta_3_ds], [0 0 0]);
c = sym_expression2value(c,[theta_1_d^2 theta_2_d^2 theta_3_d^2], [theta_1_ds theta_2_ds theta_3_ds]);
C = equationsToMatrix(c,[theta_1_ds theta_2_ds theta_3_ds]);
% adding a trajectory to the leg end




