function [Matrix] = linVelocity(position, varVector)
varVec = formula(varVector);
pos = formula(position);
vars = length(varVec);
Matrix = sym(zeros(3,vars));
temp = sym(zeros(3,1));

for i=1:vars
   for j=1:3
       temp(j,1) = functionalDerivative(pos(j),varVec(i));
   end
   Matrix(:,i)=temp;
end

