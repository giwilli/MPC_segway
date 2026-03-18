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
%%
P_Kalm = (10^6)*eye(dim_A);
Q_Kalm = 0.5*eye(dim_A);
R_Kalm = 0.5*eye(dim_C);

w_fn = sqrt(Q_Kalm)*randn(dim_A,M+1);
v_fn = sqrt(R_Kalm)*randn(dim_C,M+1);


x_real_fn = zeros(dim_A, M+1);
x_pred_fn = zeros(dim_A, M+1);
x_real_fn(:,1) = x0;
x_pred_fn(:,1) = x0;
for i = 1:M 
    y = C*x_real_fn(:,i) + v(:,i); 
    [x_pred_fn(:, i+1),P_Kalm] = Kalm_fn(sys_d.A, sys_d.B, sys_d.C, sys_d.D,x_pred_fn(:,i),P_Kalm,Q_Kalm,R_Kalm,y,-K*x_real_fn(:,i));
    x_real_fn(:,i+1) = A*x_real_fn(:,i) + B*(-K*x_real_fn(:,i)) + w(:,i);
end
hold on;
plot(t,x_pred_fn(1,:))
plot(t,x_real_fn(1,:))

legend("x_pred","x_real")