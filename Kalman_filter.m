A = sys_d.A;
B = sys_d.B;
C = sys_d.C;
D = sys_d.D;
P_Kalm = (10^6)*eye(dim_A);
Q_Kalm = 0.5*eye(dim_A);
R_Kalm = 0.5*eye(dim_C);

x_real = zeros(dim_A, M+1);
x_pred = zeros(dim_A, M+1);
x_real(:,1) = x0;
x_pred(:,1) = x0;

w = sqrt(Q_Kalm)*randn(dim_A,M+1);
v = sqrt(R_Kalm)*randn(dim_C,M+1);

for i = 1:M 
    L = (P_Kalm*C.')/(C*P_Kalm*C.' + R_Kalm);
    y = C*x_real(:,i) + v(:,i);
    x_pred(:, i+1) = (A - A*L*C)*x_pred(:,i) + (B - A*L*D)*(-K*x_real(:,i)) + A*L*y;
    P_Kalm = A*P_Kalm*A.' - A*P_Kalm*C.'/(C*P_Kalm*C.' + R_Kalm)*C*P_Kalm*A.' + Q_Kalm;
    x_real(:,i+1) = A*x_real(:,i) + B*(-K*x_real(:,i)) + w(:,i);
end


hold on;
plot(t,x_pred(1,:))
plot(t,x_real(1,:))

legend("x_pred","x_real")