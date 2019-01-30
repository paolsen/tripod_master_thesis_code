function [K] = v2zM(matrix,columns)
M = formula(matrix);
[n,m] = size(M);
Z = sym(zeros(n,columns));
for i=1:m
    Z(:,i) = M(:,i);
end
K = Z;
end


