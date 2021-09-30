% Kalman Filter
%%
%
%
%

clc;
clear;

I=1:20000;

M = csvread('r-gelio.csv');
V = csvread('V-gelio.csv');
T = csvread('t-gelio.csv');

E = eye(6);
R = [1 0 0 0 0 0;
     0 1 0 0 0 0;
     0 0 1 0 0 0;
     0 0 0 1e-10 0 0;
     0 0 0 0 1e-10 0;
     0 0 0 0 0 1e-10];

R_1 = M(1,:) * (1 + (1e8+64e5)/norm(M(1,:))) + [0 64e5 0] ;
R_2 = M(1,:) * (1 + (1e8+64e5)/norm(M(1,:))) + [0 -64e5 0] ;

State_X=[M,V];

% State_X = [T(2:size(T,1)),M(2:size(M,1),:)];
% for i=1:1:size(M,1)-1
%    State_X(i,5:7) = (M(i+1,:) - M(i,:))/(T(i+1)-T(i));  
% end

x = zeros(20000,6);
x(1,:)=ones(1,6)*1e-3;
K = zeros(6,6,20000);
K(:,:,1) = eye(6);

for i=2:20000
    dt = T(i)-T(i-1);
    x_k_1 = x(i-1,:);
    K_k_1=K(:,:,i-1);
    F = calc_F(x_k_1,dt);
    x_p = predict_x(F,x_k_1);
    K_p = predict_K(F,K_k_1);

    H = calc_H(x_k_1,R_1,R_2);
    a = calc_a(K_p,H,R);
    %y = calc_y(H);
    
    x(i,:) = refine_x(x_p,a,H,R)'; 
    K(:,:,i) = refine_K(a,H,K_p); 
    save('Kalman_filter_2.mat');
    
    
    %plot(x(1:i,1),I(1:i));
    %drawnow
    
end
plot(T,vecnorm(x(:,1:3),2,2))