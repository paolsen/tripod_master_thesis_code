% Test sytem


syms x(t) y(t) x_d(t) y_d(t) x_dd(t) y_dd(t)
syms m1 m2
syms I1 I2
syms l1 l2

B0 = DH(0,sym(pi)/2,l1,x);
B1 = DH(0,0,l2,y);

mvec = [m1,m2];

X = [x;y];
X_D = [x_d;y_d];
X_DD = [x_dd;y_dd];

cm1 = four2three(B0)*[0; 0; l1/2;1];
cm2 = four2three(B0*B1)*[0; 0; l2/2;1];


j1 = linVelocity(pos(B0),X);
j2 = linVelocity(pos(B1),X);

J_v = [j1,j2];
z_rot = sym([0;0;1]);
x_rot = sym([1;0;0]);
jw1 = v2zM(z_rot,2);
jw2 = v2zM([z_rot,z_rot],2);


J_w = [jw1,jw2];

R = [rot(GenRot(x,0,0)),rot(GenRot(y,0,0))];
I = [eye(3)*I1 eye(3)*I2];

D = InertiaMatrix(mvec,J_v,J_w,R,I);
K = 0.5*X_D.'*D*X_D;

syms g
P = g*mvec*[cm1(3,1);cm2(3,1)];

L = K-P;

Lt = timeDependant(L,X,X_D,t);


eq = DeriveEL(L,X,X_D,X_DD,t);

Eq = char(eq);
fileID = fopen('equations.txt','w');
fprintf(fileID,Eq);
fclose(fileID);
