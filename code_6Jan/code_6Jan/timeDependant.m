function [newFunc] = timeDependant(f,v,v_d,variable)
v = formula(v);
v_d = formula(v_d);
v_dt = diff(v,variable);
for i=1:length(v)
    f = subs(f,v_d(i,1),v_dt(i,1));
end
newFunc = f;
end

