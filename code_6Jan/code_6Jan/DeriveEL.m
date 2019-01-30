function [equations] = DeriveEL(Lagrangian, q, q_d, q_dd, time)
    L = formula(Lagrangian);
    q = formula(q);
    q_d = formula(q_d);
    q_dd = formula(q_dd);
    n = length(q);
    eqVector = -functionalDerivative(L,q);
    eqMatrix = functionalDerivative(L,q_d);
    q_dt = diff(q,time);
    q_ddt = diff(q_d,time);
    for i=1:n
        eqMatrix(i,1) = diff(eqMatrix(i,1),time);
        for j=1:n
            eqMatrix(i,1) = subs(eqMatrix(i,1),q_dt(j,1),q_d(j,1));
            eqMatrix(i,1) = subs(eqMatrix(i,1),q_ddt(j,1),q_dd(j,1));
        end
    end
    equations = eqMatrix + eqVector;
end

