%polar_plot_test


t = linspace(0,2*pi/6,1000);
theta = -pi/6+t;
r = 3*1./cos(theta);

%polarplot(theta,r);

x = r.*cos(theta);
y = r.*sin(theta);

axis('equal');
plot(x,y);