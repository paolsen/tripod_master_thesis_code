%two-leg-stance-EL

syms m g d(t) d(t) phi1(t) t
syms d_d(t) phi1_d(t)
syms d_dd(t) phi1_dd(t)

q = [d; phi1];
q_d = [d_d; phi1_d];
q_dd = [d_dd;phi1_dd];

P = m*g*d*cos(phi1);
K = 0.5*(m*d^2*phi1_d^2 + m*d_d^2);

L = K-P;

eq = DeriveEL(L,q,q_d,q_dd,t)