function [vector] = LegTransform(th1,th2,th3,d1,d2,d3,r)
vector = [r+(d1+d2*cos(th2)+d3*cos(th2+th3))*cos(th1);(d1+d2*cos(th2)+d3*cos(th2+th3))*sin(th1);-(d2*sin(th2)+d3*sin(th2+th3))];
end

