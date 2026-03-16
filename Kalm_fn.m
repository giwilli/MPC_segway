function [x_pred] = Kalm_fn(A_sys, B_sys, C_sys, D_sys, x_prev_real, x_prev_pred)
%UNTITLED2 Kalman filter
%   Implementation of the Kalman filter for the OTS
A = A_sys;
B = B_sys;
C = C_sys;
D = D_sys;

dim_A = size(A,1);
dim_B = size(B,2);
dim_C = size(C,1);

P_Kalm = (10^6)*eye(dim_A);
R_Kalm = 0.5*eye(dim_C);

L = (P_Kalm*C.')/(C*P_Kalm*C.' + R_Kalm);
y = C*x_prev_real;
x_pred = (A - A*L*C)*x_prev_pred + (B - A*L*D)*(-K*x_prev_real) + A*L*y;
end