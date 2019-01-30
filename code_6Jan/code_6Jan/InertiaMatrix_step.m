function [D_step] = InertiaMatrix_step(mass,LinVelJac,AngVelJac,Rotation,Inertia)
m = formula(mass);
Jv = formula(LinVelJac);
Jw = formula(AngVelJac);
R = formula(Rotation);
I = formula(Inertia);
D_step = m*Jv.'*Jv + Jw.'*R*I*R.'*Jw;
end

