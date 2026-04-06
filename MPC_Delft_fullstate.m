%% Cleanup and initialize

clc;
tbxmanager restorepath;
mpt_init;

%% Global Solver option
global MPTOPTIONS
MPTOPTIONS.modules.geometry.sets.Polyhedron.projection.method = 'mplp';
%% Create the discretized system
%clear all;
close all;
clc;
Equations;
load_Tset = false;
A = A_lin_s;
B = B_lin_s;
C = C_lin_s;
D_sys = D_lin_s;

Ts = 0.01;
sys_d = c2d(ss(A,B, C, D_sys),Ts,'zoh');

rank(ctrb(sys_d.A, sys_d.B))
rank(obsv(sys_d.A, sys_d.C))

%% Problem Fundamentals

N = 45;
dim_A = size(A,1);
dim_B = size(B,2);
dim_C = size(C,1);
Q = 1000*eye(dim_A);
R = 1*eye(dim_B);

c = [100; 7; pi/18; 100];
u_bound = 42;

%% Constraints definition and Terminal Set

model = LTISystem(sys_d);
model.x.min = -c;
model.x.max = c;
model.u.min = -u_bound;
model.u.max = u_bound;

model.x.penalty = QuadFunction(Q);
model.u.penalty = QuadFunction(R);
P = model.LQRPenalty.weight;

%% Compact MPC Formulation

T = zeros(dim_A*(N+1),dim_A); 
S = zeros(dim_A*(N+1), dim_B*N);


for i = 1:N+1
    T((i-1)*dim_A + 1 : i*dim_A, 1:dim_A) = sys_d.A^(i-1);
    if i <= N
        S(:,end-i*dim_B + 1:end -(i-1)*dim_B) = circshift(T,-dim_A*i,1)*sys_d.B;
    end
end

D = diag([1 1 1 1]);

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

if load_TSet
    Tset_Aload = load("data/Tset_A_Q1000R1.mat");
    Tset_bload = load("data/Tset_b_Q1000R1.mat");

    Tset_A = Tset_Aload.Tset_A;
    Tset_b = Tset_bload.Tset_b;
else
    Tset = model.LQRSet;
    Tset_A = Tset.A;
    Tset_b = Tset.b;

D_tilde_term = [D_terminal*sys_d.A];   %;-D_terminal*sys_d.A; zeros(1,dim_A); zeros(1,dim_A)];
E_tilde_term = [D_terminal*sys_d.B];   %;-D_terminal*sys_d.B; 0; 0];
b_tilde_term = [c_terminal];           %;c_terminal;0;0];

tmp = zeros(1,N);
tmp(:,end) = 1;
D_bar_term_temp = kron(tmp,D_tilde_term);
E_bar_term_temp = kron(tmp,E_tilde_term);
b_bar_term_temp = b_tilde_term;

%% Constraint Concatination
x0 = [10;0;0;0];

D_bar = [D_bar_temp;D_bar_term_temp];
E_bar = [E_bar_temp;E_bar_term_temp];
b_bar = [b_bar_temp;b_bar_term_temp];

G = D_bar*S_tilde + E_bar;

g = b_bar - D_bar*T_tilde*x0;

%% Closed Loop global paramters

sim_sec = 10;
t = 0:Ts:sim_sec;
M = sim_sec/Ts;
t = 0:Ts:M*Ts;
y_ref_final = 1;
y_constant = ones(1,M+1)* y_ref_final;
y_square = square (pi*t/10);
y_linear = linspace(0,y_ref_final,M);
y_ref = y_constant;%[linspace(0,y_ref_final,M)];
%% Closed Loop MPC

x_mpc = zeros(dim_A, (M+1));
x_mpc(:,1) = x0;
u_mpc_log = zeros(dim_B, M+1);
u_mpc_log(:,1) = 0;

for i = 1:M
    disp(i);
    x0 = x_mpc(:,i);
    h = S.'*Q_bar*T*x0;
    g = b_bar - D_bar*T_tilde*x0;
    u = quadprog(H,h, G,g);
    u_mpc_log(:,i) = u(1,:);
    x_mpc(:,i+1) = sys_d.A*x_mpc(:,i) + sys_d.B*u(1);
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

% subplot(2,1,1);
% plot(t,x_mpc(1,:))
% title('State Trajectories (x)')
% legend('x1')
% subplot(2,1,2);
% plot(t,u_mpc_log(1,:))
% %plot(t, x_mpc(1,:))
% %plot(t, x_mpc(1,:), t, x_lqr);
% title('Control Input (u)');
% legend('u');
%% Plot the different horizons
% subplot(2,1,1)
% hold on
% plot(t,x_N1(1,:))
% plot(t,x_N5(1,:))
% plot(t,x_N10(1,:))
% plot(t,x_N25(1,:))
% plot(t,x_N50(1,:))
% plot(t,x_N75(1,:))
% plot(t,x_N100(1,:))
% legend('1', '5', '10', '25', '50', '75', '100');
% title('State Trajectories (x)');

