function [x_ref,u_ref] = OTS(y_ref, H_sel, A_sys,B_sys,C_sys,Bd_sys, Cd_sys,d_hat,D,c,u_bound,max_it,H_ots,h_ots)
%OTS OTS function
%   Function to calculate the optimal target selection
dim_A = size(A_sys,1);
dim_B = size(B_sys,2);
dim_C = size(C_sys,1);
dim_HC = size(H_sel*C_sys,1);
A_aug = [(eye(dim_A)-A_sys) -B_sys;
        H_sel*C_sys zeros(dim_HC,dim_B)];
HC_d = H_sel*Cd_sys;
b_aug = [Bd_sys*d_hat;(y_ref - HC_d*d_hat)];

f = [c; u_bound];
D_aug = [D; zeros(dim_B, dim_A)];
E_aug = [zeros(dim_A,dim_B); ones(dim_B)];
F_aug = [D_aug, E_aug;-D_aug, -E_aug];
f_aug = [f; f];

options = optimoptions('quadprog', 'MaxIterations', max_it);
lb = [];
ub = [];
x0 = [];
state_aug = quadprog(H_ots, h_ots, F_aug, f_aug, A_aug, b_aug,lb,ub,x0, options);

x_ref = state_aug(1:4,1);
u_ref = state_aug(5:end,1);
end