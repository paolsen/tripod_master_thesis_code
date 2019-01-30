%test_animated plot:

%th = linspace(0,pi,1000);
%link_1_end = plot(sin(th(1)), 2-cos(th(1)),'o','MarkerFaceColor','red');
%link_1_end = animatedline('Color','blue');
link_2_end = animatedline('Color','red');
set(gca,'xLim',[-3, 3],'yLim',[-3,3]);
axis('equal');
title('double_pendulum');
legend('link_1','link_2');
grid on;
hold('on')
for i=1:length(T)-1
   plot(sin(th_1(i)),2-cos(th_1(i)));
   addpoints(link_2_end, sin(th_1(i))+sin(th_1(i)+th_2(i)),2-cos(th_1(i))-cos(th_1(i)+th_2(i)));
   drawnow
   pause(0.001);% comment out for full speed
end
hold('off');