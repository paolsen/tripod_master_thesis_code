function [matrix3] = rot(fourXfour)
M = formula(fourXfour);
matrix3 = M([1,2,3],[1,2,3]);
end

