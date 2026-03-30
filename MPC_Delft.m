%% Cleanup and initialize
% 
% clc;
%tbxmanager restorepath;
%mpt_init;
% 
% %% Global Solver option
% global mptOptions
% mptOptions.verbose = 1;

%% Create the discretized system
clear all;
close all;
clc;
Equations;
A = A_lin_s;
B = B_lin_s;
C = C_lin_s;
D_sys = D_lin_s;

Ts = 0.01;
sys_d = c2d(ss(A,B, C, D_sys),Ts,'zoh');

sys_d.A
sys_d.B
sys_d.C
sys_d.D

%rank(ctrb(sys_d.A, sys_d.B))
%rank(obsv(sys_d.A, sys_d.C))

%% Problem Fundamentals

N = 100;
dim_A = size(A,1);
dim_B = size(B,2);
dim_C = size(C,1);
Q = 1000*eye(dim_A);
R = 1*eye(dim_B);

c = [Inf; 7; pi/18; Inf];
u_bound = 42;

%% Constraints definition and Terminal Set

% model = LTISystem(sys_d);
% model.x.min = -c;
% model.x.max = c;
% model.u.min = -u_bound;
% model.u.max = u_bound;
% 
% model.x.penalty = QuadFunction(Q);
% model.u.penalty = QuadFunction(R);
% P = model.LQRPenalty.weight;

%% Compact MPC Formulation

T = zeros(dim_A*(N+1),dim_A); 
S = zeros(dim_A*(N+1), dim_B*N);


for i = 1:N+1
    T((i-1)*dim_A + 1 : i*dim_A, 1:dim_A) = sys_d.A^(i-1);
    if i <= N
        S(:,end-i*dim_B + 1:end -(i-1)*dim_B) = circshift(T,-dim_A*i,1)*sys_d.B;
    end
end

D = diag([0 1 1 0]);

D_tilde = [D*sys_d.A;-D*sys_d.A; zeros(1,dim_A); zeros(1,dim_A)];
E_tilde = [D*sys_d.B;-D*sys_d.B; 1; -1];
b_tilde =[c;c;u_bound;u_bound];
D_bar_temp = kron(eye(N),D_tilde);
E_bar_temp = kron(eye(N),E_tilde);
b_bar_temp = repmat(b_tilde,[N,1]);
T_tilde = T(1:end-dim_A,:);
S_tilde = S(1:end-dim_A,:);

%% Objective function weights (Compact Form)


P = load("P.mat");
P = P.P;

Q_bar = blkdiag(kron(eye(N),Q), P);
R_bar = kron(eye(N),R);

H = S.'*Q_bar*S + R_bar;
H = (H+H')/2;

%h = S.'*Q_bar*T*x0;
%% Terminal Constraint Formulation

%Tset = model.LQRSet;
%D_terminal = Tset.A;
%c_terminal = Tset.b;

Tset_A = load("TS_A.mat");
Tset_b = load("TS_b.mat");

D_terminal = Tset_A.TS_A;
c_terminal = Tset_b.TS_b;

D_tilde_term = [D_terminal*sys_d.A];   %;-D_terminal*sys_d.A; zeros(1,dim_A); zeros(1,dim_A)];
E_tilde_term = [D_terminal*sys_d.B];   %;-D_terminal*sys_d.B; 0; 0];
b_tilde_term = [c_terminal];           %;c_terminal;0;0];

tmp = zeros(1,N);
tmp(:,end) = 1;
D_bar_term_temp = kron(tmp,D_tilde_term);
E_bar_term_temp = kron(tmp,E_tilde_term);
b_bar_term_temp = b_tilde_term;

%% Constraint Concatination
x0 = [0;0;0;0];

D_bar = [D_bar_temp;D_bar_term_temp];
E_bar = [E_bar_temp;E_bar_term_temp];
b_bar = [b_bar_temp;b_bar_term_temp];

G = D_bar*S_tilde + E_bar;

g = b_bar - D_bar*T_tilde*x0;

%% Closed Loop global paramters

M = 600;
t = 0:Ts:M*Ts;
y_ref_final = 1;
y_constant = ones(1,M+1)* y_ref_final;
y_square = square (pi*t/10);
y_linear = linspace(0,y_ref_final,M);
y_ref = y_constant;%[linspace(0,y_ref_final,M)];
%% Closed Loop MPC

x_mpc = zeros(dim_A, (M+1));
y_mpc = zeros(dim_C, (M+1));
x_mpc(:,1) = x0;
y_mpc(:,1) = 0;
u_mpc_log = zeros(dim_B, M+1);
u_mpc_log(:,1) = 0;
n_d = 1;

kalman_log = zeros(dim_A+n_d, (M+1));
x_ref_normlog = zeros(dim_A, (M+1));

d = 0.01;

q_std = 0.00001;
pos_std = 0.002;
vel_std = 0.002;
angv_std = 0.02;

P_Kalm =  (1e-3)*eye(dim_A+n_d);
Q_Kalm = q_std^2*eye(dim_A+n_d);
R_Kalm = diag([pos_std^2 vel_std^2 angv_std^2]);

w = sqrt(Q_Kalm)*randn(dim_A+n_d,M+1);
v = sqrt(R_Kalm)*randn(dim_C,M+1);

%plot(w)

x_pred = [zeros(dim_A,1); zeros(n_d)];

H_sel = [1 0 0];
B_d = [0;0;0;0];
C_d_sys = [0;0;1];
H_aug = diag([0,0,0,0,1]);
h_aug = zeros(5,1);

%y = sys_d.C * x0 + C_d_sys*d + v(:,1);

A_kalm = [sys_d.A B_d; zeros(n_d,dim_A) eye(n_d)];
B_kalm = [sys_d.B; zeros(n_d)];
C_kalm = [sys_d.C, C_d_sys];


%[x_pred, P_Kalm] = Kalm_fn(A_kalm, B_kalm, C_kalm, sys_d.D, x_pred,P_Kalm,Q_Kalm,R_Kalm,y,u_mpc_log(:,1));
%d_hat = x_pred(end-n_d+1:end);

for i = 1:M
    disp(i)
    x0 = x_pred(1:4);
    d_hat = x_pred(end,1);

    % OTS
    [x_ref, u_ref] = OTS(y_ref(i),H_sel,sys_d.A,sys_d.B,sys_d.C,B_d,C_d_sys,d_hat,D,c,u_bound,200,H_aug,h_aug);
    x_ref_bar = repmat(x_ref,[N+1,1]);
    u_ref_bar = repmat(u_ref,[N,1]);
    x_ref_normlog(:,i) = norm(x0 - x_ref);
    
    % Update of the terminal constraint
    c_terminal = Tset_b.TS_b + D_terminal*x_ref;
    b_tilde_term = c_terminal;
    b_bar_term_temp = b_tilde_term;
    b_bar = [b_bar_temp;b_bar_term_temp];
    g = b_bar - D_bar*T_tilde*x0;
    
    % Update of the cost function
    h = S.'*Q_bar*(T*x0-x_ref_bar) - R_bar*u_ref_bar;
    u = quadprog(H,h, G,g);
    u_mpc_log(:,i) = u(1,:);

    x_mpc(:,i+1) = sys_d.A*x_mpc(:,i) + sys_d.B*u(1) + B_d*d + w(1:end-n_d,i);
    y = sys_d.C * x_mpc(:,i) + C_d_sys*d+ v(:,i);
    y_mpc(:,i) = y;
    
    [x_pred, P_Kalm] = Kalm_fn(A_kalm, B_kalm, C_kalm, sys_d.D,x_pred,P_Kalm,Q_Kalm,R_Kalm,y,u(1));
    %x_pred = [x_mpc(:,i+1);d];
    kalman_log(:,i) = x_pred;
end
%% LQR for reference

% [K,S,P] = lqr(sys_d,Q,R);
% x_lqr = zeros(dim_A, (M+1));
% x_lqr(:,1) = x0;
% u_lqr_log = zeros(dim_B, M+1);
% u_lqr_log(:,1) = 0;
% 
% for i = 1:M
%     %disp(i);
%     u = -K*x_lqr(:,i);
%     x_lqr(:,i+1) = (sys_d.A*x_lqr(:,i) + sys_d.B*u);
%     u_lqr_log(:,i) = u;
% end
%% Plotting the MPC and LQR Response

subplot(4,1,1); % Top plot
plot(t, y_mpc(1,:))
subplot(4,1,2);
hold on
plot(t,kalman_log(1,:))
plot(t,x_mpc(1,:))
subplot(4,1,3);
plot(t,kalman_log(5,:))
subplot(4,1,4);
plot(t,y_ref)
%plot(t, x_mpc(1,:))
%plot(t, x_mpc(1,:), t, x_lqr);
title('State Trajectories (x)');
legend('Reference');

%%
% subplot(2,1,2); % Bottom plot
% plot(t, u_mpc_log, t, u_lqr_log);
% title('Control Inputs (u)');
% legend('MPC', 'LQR');
% 
% %% OTS
% H_sel = [1 0];
% dim_HC = size(H_sel*sys_d.C,1);
% A_aug = [(eye(dim_A)-sys_d.A) -sys_d.B;
%         H_sel*sys_d.C zeros(dim_HC,dim_B)];
% %det(A_aug);
% d_hat = 0;
% y_ref = 1;
% B_d = [0;0;0;0];
% C_d_sys = [0;0.1];
% HC_d = H_sel*C_d_sys;
% b_aug = [B_d*d_hat;(y_ref - HC_d*d_hat)];
% 
% % sys_d.C;
% % A_distrej = [(eye(dim_A)-sys_d.A) -sys_d.B;
% %         H_sel*sys_d.C C_d];
% % 
% % rank(A_distrej)
% % rank(A_aug)
% 
% D = diag([1 1 1 1]);
% c = [1000; 7; pi/18; 1000];
% u_bound = 42;
% f = [c; u_bound];
% D_aug = [D; zeros(dim_B, dim_A)];
% E_aug = [zeros(dim_A,dim_B); ones(dim_B)];
% F_aug = [D_aug, E_aug;-D_aug, -E_aug];
% f_aug = [f; f];
% 
% H_aug = diag([0,0,0,0,1]);%eye(5);%
% h_aug = zeros(5,1);
% 
% options = optimoptions('quadprog', 'MaxIterations', 200);
% lb = [];
% ub = [];
% x0 = [];
% state_aug = quadprog(H_aug, h_aug, F_aug, f_aug, A_aug, b_aug,lb,ub,x0, options);
% 
% x_ref = state_aug(1:4,1);
% u_ref = state_aug(5:end,1);
% 
% %%
% [x_ref_fn,u_ref_fn] = OTS(y_ref,H_sel,sys_d.A,sys_d.B,sys_d.C,B_d,C_d_sys,d_hat,D,c,u_ref,200,H_aug,h_aug)
% % Define the constraints for x_ref and u_ref
% % Define the constraints for y_ref
% % Solve the Optimal Control Problem
