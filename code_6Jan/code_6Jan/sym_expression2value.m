function [value] = sym_expression2value(expression, sym_vector, number_vector)
for i=1:length(sym_vector)
    expression = subs(expression,sym_vector(i),number_vector(i));
end
value = expression;
end

