% simple model pendulum motion
steps = 10000;
t = linspace(0,5,steps);
d = zeros(1,steps);
d_d = zeros(1,steps);
d_dd = zeros(1,steps);
phi = zeros(1,steps);
phi_d = zeros(1,steps);
phi_dd = zeros(1,steps);
f_perp = 8.923035;
f_list = [8.7 8.8 8.9 8.95];
f_d = 1;
t_1 = 0.5;
d(1) = 0.9;
phi(1) = pi/4;
m = 1;
g = 9.81;
h = (t(steps)-t(1))/(length(t))


for j = 1:4
    for i=1:(length(t)-1)
        if(t(i) < t_1)
            %f = f_perp;
            f = f_list(j);
        else
            f = 0;
        end
    
        phi_dd(i) = (-2*m*g*d(i)*d_d(i)*phi_d(i) + g*m*d(i)*sin(phi(i))-f*d(i))/(m*d(i)^2);
        phi_d(i+1) = phi_d(i) + phi_dd(i)*h;
        phi(i+1) = phi(i) + phi_d(i)*h;
        %d_dd(i) = -g*cos(phi(i))+d(i)*phi(i) + d(i)*phi_d(i)^2 + (f_d*m*cos(phi(i))*g)/m;
        %d_d(i+1) = d_d(i) + d_dd(i)*h;
        %d(i+1) = d(i) + d_d(i)*h;
        d(i+1)= d(i);
    end
    %figure(1)
    %plot(t,phi);
    plot(phi,phi_d,'LineWidth',2);
    hold('on')
end


figure(1)
%plot(t,phi);
%hold('on');
%plot(t,0.25*pi*ones(1,length(t)),'color',[0,0,0]+0.7);
%plot(t,-0.25*pi*ones(1,length(t)),'color',[0,0,0]+0.7);
%plot(t,zeros(1,length(t)),'color',[0,0,0]);
%yticks([-pi -1/2*pi -1/4*pi 0 1/4*pi 1/2*pi pi]);
%yticklabels({'-\pi','-1/2\pi','-1/4\pi','0','1/4\pi','1/2\pi','\pi'});
%axis([0 5 -pi pi]);
legend('','f_p = 8.7','f_p = 8.8','f_p = 8.9', 'f_p = 8.95')
%xlabel('time (seconds)')
%ylabel('Angle: \phi_1')
%figure(2)
%plot(t,phi_d);
%xlabel('time')
%yl = ylabel('Angular velocity: $\dot{\phi}_1$');
%set(yl, 'Interpreter', 'latex');
%hl = legend('$\dot{\phi}_1$');
%set(hl, 'Interpreter', 'latex');
% multiple plots
    
