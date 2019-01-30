function [D] = InertiaMatrix(massVector,LinVelJac,AngVelJac,Rotations,InertiaM)
m = formula(massVector);
Jv = formula(LinVelJac);
Jw = formula(AngVelJac);
R = formula(Rotations);
I = formula(InertiaM);
n = length(massVector);
it = 0;
it2 = 0;
E = sym(zeros(n));
%m(1)
%Jv(:,ones(1,n)*it + 1:n)
%Jw(:,ones(1,n)*it+ 1:n)
%R(:,ones(1,3)*it2+[1,2,3])
%I(:,ones(1,3)*it2+[1,2,3])
for link=1:n
    E = E + InertiaMatrix_step(m(link),Jv(:,ones(1,n)*it + [1:n]),Jw(:,ones(1,n)*it+ [1:n]),R(:,ones(1,3)*it2+[1,2,3]),I(:,ones(1,3)*it2+[1,2,3]));
    it = it+n;
    it2= it2+3;
end
D = E;
end

