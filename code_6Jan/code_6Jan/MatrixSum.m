function [Matrix] = MatrixSum(nXnm)
M = formula(nXnm);
[n,nm] = size(M);
m = nm/n;
S = superEye(n,m);
Matrix = M*S;
end

