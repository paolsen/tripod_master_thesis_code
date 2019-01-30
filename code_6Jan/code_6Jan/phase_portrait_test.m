%phase plane test

g = 9.81;
m = 1;
d = 0.9;
%f = @(t,Y) [Y(2); -sin(Y(1))];
f = @(t,Y) [Y(2); g/d*sin(Y(1))];
%y1 = linspace(-2,8,20);
%y2 = linspace(-2,2,20);
y1 = linspace(-2*pi,2*pi,20);
y2 = linspace(-6,6,20);

% creates two matrices one for all the x-values on the grid, and one for
% all the y-values on the grid. Note that x and y are matrices of the same
% size and shape, in this case 20 rows and 20 columns
[x,y] = meshgrid(y1,y2);
u = zeros(size(x));
v = zeros(size(x));

% we can use a single loop over each element to compute the derivatives at
% each point (y1, y2)
t=0; % we want the derivatives at each point at t=0, i.e. the starting time
for i = 1:numel(x)
    Yprime = f(t,[x(i); y(i)]);
    u(i) = Yprime(1);
    v(i) = Yprime(2);
end
for i = 1:numel(x)
    Vmod = sqrt(u(i)^2 + v(i)^2);
    u(i) = u(i)/Vmod;
    v(i) = v(i)/Vmod;
end
quiver(x,y,u,v,'color',[0,0,0]+0.7); figure(gcf)
xticks([-pi  0  pi]);
xticklabels({'-\pi','0','\pi'});
xlabel('Angle: \phi_1')
ylab = ylabel('Angular velocity: $\dot{\phi}_1$');
set(ylab, 'Interpreter', 'latex');
axis tight equal;
axis([-2*pi, 2*pi, -6, 6]);

hold('on');
%for y20 = [0 0.5 1 1.5 2 2.5]
%    [ts,ys] = ode45(f,[0,50],[0;y20]);
%    plot(ys(:,1),ys(:,2),'b')
    %plot(ys(1,1),ys(1,2),'bo') % starting point
    %plot(ys(end,1),ys(end,2),'ks') % ending point
%end
%hold off
