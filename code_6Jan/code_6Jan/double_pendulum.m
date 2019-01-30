%test feil i ligninnger
%syms g L_1 L_2 m_1 m_2 theta_1(t) theta_2(t) theta_1_d(t) theta_2_d(t) theta_1_dd(t) theta_2_dd(t) t
syms theta_1(t) theta_2(t) theta_1_d(t) theta_2_d(t) theta_1_dd(t) theta_2_dd(t) t
g = 10;
L_1 = 1;
L_2 = 1;
m_1 = 1;
m_2 = 2;

K = 0.25 * theta_1_d(t)^2*m_1*L_1*(cos(theta_1)+sin(theta_1)) + (0.5*theta_1_d^2*L_1*(cos(theta_1)+sin(theta_1))+0.25*theta_2_d(t)^2*L_2*(cos(theta_2)+sin(theta_2)))*m_2;
P = m_1*g*(L_1+L_2-0.5*L_1*cos(theta_1))+m_2*g*(L_1+L_2-L_1*cos(theta_1)-0.5*L_2*cos(theta_1 + theta_2));

L = K -P;

q = [theta_1; theta_2];
q_d = [theta_1_d; theta_2_d];
q_dd = [theta_1_dd; theta_2_dd];

eqs = DeriveEL(L,q,q_d,q_dd,t);

eq_1 = eqs(1);
alpha_1 = subs(subs((eq_1- subs(eq_1,theta_1_dd,0)),theta_1_dd,1),theta_1_d,0);
beta_1 = subs(subs((eq_1 - subs(eq_1,theta_1_d,0)),theta_1_d,1),theta_1_dd,0);
gamma_1 = subs(subs(subs(eq_1,theta_1_dd,0),theta_1_d,0),theta_2_d,0);
eq_2 = eqs(2); 
alpha_2 = subs(subs((eq_2 - subs(eq_2,theta_2_dd,0)),theta_2_dd,1),theta_2_d,0);
beta_2 = subs(subs((eq_2 - subs(eq_2,theta_2_d,0)),theta_2_d,1),theta_2_dd,0);
gamma_2 = subs(subs(subs(eq_2,theta_2_dd,0),theta_2_d,0),theta_1_d,0);



Alpha = [alpha_1;alpha_2];
Beta = [beta_1;beta_2];
Gamma = [gamma_1;gamma_2];

T = linspace(0,10,1000);
th_1 = zeros(1,1000);
th_1(1) = pi/2;
th_2 = zeros(1,1000);
th_2(1) = 0;
th_1_d = zeros(1,1000);
th_1_d(1) = 0;
th_2_d = zeros(1,1000);
th_2_d(1) = 0;
th_1_dd = zeros(1,1000);
th_2_dd = zeros(1,1000);
a_1 = zeros(1,1000);
a_2 = zeros(1,1000);
b_1 = zeros(1,1000);
b_2 = zeros(1,1000);
g_1 = zeros(1,1000);
g_2 = zeros(1,1000);
dt = 0.002;

link_1_end = plot(sin(th_1(1)),2-cos(th_1(1)),'o','MarkerFaceColor','red');
%link_1_end = animatedline('Color','blue');
link_2_end = animatedline('Color','red');
set(gca,'xLim',[-10, 10],'yLim',[-10,10]);
axis('equal');
title('double_pendulum');
legend('link_1','link_2');
grid on;
hold('on')
for i=1:length(T)-1
   a_1(i) = subs(subs(alpha_1,theta_1,th_1(i)),theta_2,th_2(i));
   a_2(i) = subs(subs(alpha_2,theta_1,th_1(i)),theta_2,th_2(i));
   b_1(i) = subs(subs(beta_1,theta_1,th_1(i)),theta_2,th_2(i));
   b_2(i) = subs(subs(beta_2,theta_1,th_1(i)),theta_2,th_2(i));
   g_1(i) = subs(subs(gamma_1,theta_1,th_1(i)),theta_2,th_2(i));
   g_2(i) = subs(subs(gamma_2,theta_1,th_1(i)),theta_2,th_2(i));
   th_1_dd(i+1) = (b_1*th_1_d(i)^2 + g_1)/a_1;
   th_2_dd(i+1) = (b_2*th_2_d(i)^2 + g_2)/a_2;
   th_1_d(i+1) =th_1_d(i) + th_1_dd(i+1)*dt;
   th_2_d(i+1) =th_2_d(i) + th_2_dd(i+1)*dt;
   th_1(i+1) =th_1(i) + th_1_d(i+1)*dt;
   th_2(i+1) =th_2(i) + th_2_d(i+1)*dt;
end
for i=1:length(T)-1
   link_1_end.XData = sin(th_1(i));
   link_1_end.YData = 2-cos(th_1(i));
   addpoints(link_2_end, sin(th_1(i))+sin(th_1(i)+th_2(i)),2-cos(th_1(i))-cos(th_1(i)+th_2(i)));
   drawnow
   %pause(0.001);% comment out for full speed
end
hold('off');