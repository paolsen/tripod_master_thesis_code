% simulation of moving platform while 3-leg stance
% four torques

%coordinate system pointing x = -c, z = -g

%configuration polygon
A_pos = [0.5*sqrt(1-0.25); -0.5;0];
B_pos = [0.5*sqrt(1-0.25); 0.5;0];
C_pos = [-0.5*sqrt(1-0.25);0;0];

CM_x = zeros(1,100);
CM_y = zeros(1,100);
CM_z = zeros(1,100);

m = 1;
g = 10;
f = 1;

forces = zeros(3,100);
ext_force = zeros(3,100);
a_forces = zeros(3,100);
b_forces = zeros(3,100);
c_forces = zeros(3,100);

t = linspace(0,10,100);
for i=1:99
    CM_x(i+1) = CM_x(i)+0.006;
    %CM_y(i+1) = CM_y(i)+0.001;
    [A,B,C] = CM_pos_to_normal_force([CM_x(i);CM_y(i)],A_pos,B_pos,C_pos);
    forces(1,i) = A*m*g;
    forces(2,i) = B*m*g;
    forces(3,i) = C*m*g;
    
end

plot(t,forces(1,:),t,forces(2,:),t,forces(3,:));
legend(['A'],['B'],['C']);
title('leg forces');

%plot3(forces(1,:),forces(2,:),forces(3,:));

