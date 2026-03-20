function [x_pred,P_Kalm] = Kalm_fn(A_sys, B_sys, C_sys, D_sys, x_prev_pred,P,Q,R,y_prev,u_prev)
%KALMAN_FILTER
%   Implementation of the Kalman filter for the OTS
A = A_sys;
B = B_sys;
C = C_sys;
D = D_sys;

P_Kalm = P;
Q_Kalm = Q;
R_Kalm = R;
y = y_prev;
u = u_prev;

L = (P_Kalm*C.')/(C*P_Kalm*C.' + R_Kalm);
x_pred = (A - A*L*C)*x_prev_pred + (B - A*L*D)*u + A*L*y;
P_Kalm = A*P_Kalm*A.' - A*P_Kalm*C.'/(C*P_Kalm*C.' + R_Kalm)*C*P_Kalm*A.' + Q_Kalm;
end