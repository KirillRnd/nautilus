function [H] = calc_H(State,R_1,R_2)
r=State(1:3);
V=State(4:6);
r_12 = r - R_1;
r_23 = r - R_2;

v_12 = r_12/norm(r_12)*dot(V,r_12/norm(r_12));
v_23 = r_23/norm(r_23)*dot(V,r_23/norm(r_23));

%Производная радиальной скорости
dD = 1/2*[v_12*r_12'/norm(r_12)+v_23*r_23'/norm(r_23)];

H = 1/2*[v_12'/norm(r_12)+v_23'/norm(r_23)-dD*(r_12'/norm(r_12)+r_23'/norm(r_23));
         r_12'/norm(r_12)+r_23'/norm(r_23)];

H = diag(H);     
     
     
end