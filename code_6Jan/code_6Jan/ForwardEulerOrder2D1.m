function [q,q_d,q_dd] = ForwardEulerOrder2D1(acceleration_eq,Q,Q_dot,ext_force_sym,ext_force,parameters_sym,parameters_float,x_0,v_0,t_0,t_1,steps)
%acceleration is a symbolic function
%Q is symbolic q
%Q_dot is symbolic q_dot
%ext_force_sym is the symbolic external force
%ext_force is the array of force-values as long as steps
%parameters_sym are the array of symbolic parameters
%parameters float are the array of float values coresponding to the latter
%x_0 is the initial position of q, float
%v_0 is the initial velocity (q_dot0), float
%t_0 is the start time, float
%t_1 is the stop time, float
%steps is the number of elements, int 
q = zeros(1,steps);
q(1) = x_0;
q_d = zeros(1,steps);
q_d(1) = v_0;
q_dd = zeros(1,steps);
h = (t_1-t_0)/steps;
acceleration_eq = sym_expression2value(acceleration_eq, parameters_sym,parameters_float);

for i=1:steps-1
    q_dd(i) = sym_expression2value(acceleration_eq, [Q; Q_dot; ext_force_sym],[q(i); q_d(i); ext_force(i)]);
    q_d(i+1) = q_d(i) + q_dd(i)*h;
    q(i+1)= q(i) + q_d(i)*h;
end

