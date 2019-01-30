function [distance] = dist_line_point(L1,L2,P)
distance = abs((L2(2)-L1(2))*P(1)-(L2(1)-L1(1))*P(2) +L2(1)*L1(2)-L2(2)*L1(1))/sqrt((L2(2)-L1(2))^2+(L2(1)-L1(1))^2);
end