function [n_A,n_B,n_C] = CM_pos_to_normal_force(CM,A_pos,B_pos,C_pos)
h_a = dist_line_point(C_pos,B_pos,A_pos);
h_b = dist_line_point(C_pos,A_pos,B_pos);
h_c = dist_line_point(A_pos,B_pos,C_pos);
d_a = dist_line_point(C_pos,B_pos,CM);
d_b = dist_line_point(C_pos,A_pos,CM);
d_c = dist_line_point(A_pos,B_pos,CM);
if(d_a < h_a)
    n_A = d_a/h_a;
else
    n_A = 0;
end

if(d_b < h_b)
    n_B = d_b/h_b;
else
    n_B = 0;
end

if(d_c < h_c)
    n_C = d_c/h_c;
else
    n_C = 0;
end
end

