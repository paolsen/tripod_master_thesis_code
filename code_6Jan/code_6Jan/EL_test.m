%euler-lagrange test script:


syms a1(t) a2(t) a1_d(t) a2_d(t) a1_dd(t) a2_dd(t) t

q = [a1; a2];
q_d = [a1_d; a2_d];
q_dd = [a1_dd; a2_dd];

K = q_d.'*0.5*eye(2)*q_d;
P = [1 2] * q;

L = K-P;
%eq = DeriveEL(L,q,q_d, q_dd, t);

syms b1 b2
Eq = [2*b1 + 3*b2; 4*b1];

M = equationsToMatrix(Eq, [b1 b2])

