function [I] = superEye(identityDim, copies)
colLength = copies*identityDim;
I = eye(identityDim);
for i=1:(colLength*identityDim)
    I(1+mod(i,colLength),1+mod(i,identityDim)) = 1;
end

