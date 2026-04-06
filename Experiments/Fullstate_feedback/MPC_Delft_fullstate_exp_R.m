%% Cleanup and initialize

clc;
tbxmanager restorepath;
mpt_init;

%% Global Solver option
global MPTOPTIONS
MPTOPTIONS.modules.geometry.sets.Polyhedron.projection.method = 'mplp';
% global mptOptions
% mptOptions.verbose = 1;
% mptOptions.infbound = 500; % Increase from the default 100
%% Create the discretized system
clear all;
close all;
clc;
run('../../Equations.m');
load_TSet = false;
A = A_lin_s;
B = B_lin_s;
C = C_lin_s;
D_sys = D_lin_s;

Ts = 0.01;
sys_d = c2d(ss(A,B, C, D_sys),Ts,'zoh');

rank(ctrb(sys_d.A, sys_d.B));
rank(obsv(sys_d.A, sys_d.C));

dim_A = size(A,1);
dim_B = size(B,2);
dim_C = size(C,1);


c = [100; 7; pi/18; 100];
u_bound = 42;

model = LTISystem(sys_d);
model.x.min = -c;
model.x.max = c;
model.u.min = -u_bound;
model.u.max = u_bound;


%P = P*2;


%% Running for different N
R_values = [0.000001,0.00001,0.0001,0.001, 0.01];
N = 45;

% Initialize a structure to hold your results
results = struct();

for j = 1:length(R_values)
    R_val = R_values(j);
    disp('Current R is: ');
    disp(R_val);
    

    tStart = tic;


    Q = 1*eye(dim_A);
    R = R_val*eye(dim_B);
    model.x.penalty = QuadFunction(Q);
    model.u.penalty = QuadFunction(R);
    %P = model.LQRPenalty.weight;
    [P,K,L] = idare(sys_d.A,sys_d.B,Q,R);

    if load_TSet
        Tset_Aload = load("../../data/Tset_A_Q1000R1.mat");
        Tset_bload = load("../../data/Tset_A_Q1000R1.mat");
    
        Tset_A = Tset_Aload.Tset_A_new;
        Tset_b = Tset_bload.Tset_b_new;
    else
        Tset = model.LQRSet;
        Tset_A = Tset.A;
        Tset_b = Tset.b;
        
    end
    D_terminal = Tset_A;
    c_terminal = Tset_b;

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
    
    
    %P = load("P.mat");
    %P = P.P;
    
    Q_bar = blkdiag(kron(eye(N),Q), P);
    R_bar = kron(eye(N),R);
    
    H = S.'*Q_bar*S + R_bar;
    H = (H+H')/2;
    
    %h = S.'*Q_bar*T*x0;
   
    
    D_tilde_term = [D_terminal*sys_d.A];   %;-D_terminal*sys_d.A; zeros(1,dim_A); zeros(1,dim_A)];
    E_tilde_term = [D_terminal*sys_d.B];   %;-D_terminal*sys_d.B; 0; 0];
    b_tilde_term = [c_terminal];           %;c_terminal;0;0];
    
    tmp = zeros(1,N);
    tmp(:,end) = 1;
    D_bar_term_temp = kron(tmp,D_tilde_term);
    E_bar_term_temp = kron(tmp,E_tilde_term);
    b_bar_term_temp = b_tilde_term;
    
    x0 = [1;0;0;0];
    
    D_bar = [D_bar_temp;D_bar_term_temp];
    E_bar = [E_bar_temp;E_bar_term_temp];
    b_bar = [b_bar_temp;b_bar_term_temp];
    
    G = D_bar*S_tilde + E_bar;
    
    g = b_bar - D_bar*T_tilde*x0;
    
    
    
    sim_sec = 10;
    t = 0:Ts:sim_sec;
    M = sim_sec/Ts;
    t = 0:Ts:M*Ts;
    y_ref_final = 1;
    y_constant = ones(1,M+1)* y_ref_final;
    y_square = square(pi*t/10);
    y_sine = sin(t);
    y_linear = linspace(0,y_ref_final,M);
    y_ref = y_constant;%[linspace(0,y_ref_final,M)];

    x_mpc = zeros(dim_A, (M+1));
    x_mpc(:,1) = x0;
    u_mpc_log = zeros(dim_B, M+1);
    u_mpc_log(:,1) = 0;
    SC_log = zeros(1,M);
    TC_log = zeros(1,M);
    
    for i = 1:M
        disp(i);
        x0 = x_mpc(:,i);
        h = S.'*Q_bar*T*x0;
        g = b_bar - D_bar*T_tilde*x0;
        u = quadprog(H,h, G,g);
        u_mpc_log(:,i) = u(1,:);
        x_mpc(:,i+1) = sys_d.A*x_mpc(:,i) + sys_d.B*u(1);
        % --- COST LOGGING ---
        % Reconstruct the full predicted trajectory: X_pred = [x_0; x_1; ...; x_N]
        X_pred = T*x0 + S*u;
        U_pred = u;
        
        % 1. TERMINAL COST (TC)
        % Extract the error at the final prediction step, x_N
        x_N_pred = X_pred(end-dim_A+1 : end);
        TC_log(i) = x_N_pred.' * P * x_N_pred;
        
        X_pred_stages = X_pred(dim_A+1 : end-dim_A); 
        
        % Create the Q matrix for just the N-1 intermediate stages
        Q_stages = kron(eye(N-1), Q);
        
        % Calculate cumulative state and input costs over the horizon
        SC_states = X_pred_stages.' * Q_stages * X_pred_stages;
        SC_inputs = U_pred.' * R_bar * U_pred;
       
        SC_log(i) = SC_states + SC_inputs;
        % --------------------
    end 
    time_passed = toc(tStart);
    % Create the dynamic field name (e.g., 'x_10')
    expName = sprintf('exp_%d', j); 
    
    % Save x and u inside that specific experiment's field
    results.(expName).x = x_mpc;
    results.(expName).u = u_mpc_log;
    results.(expName).SC = SC_log;
    results.(expName).TC = TC_log;
    results.(expName).time = time_passed;
end
save('results_R.mat', 'results');
%% Plot the different horizons
% raw_fields = fieldnames(results); 
% dynamic_labels = strrep(raw_fields, 'exp_', '');
% ss_bounds = ones(1,M+1) * 0.02;
% max_overshoot = ones(1,M+1)* 0.1;
% 
% subplot(1,1,1)
% hold on
% grid on
% structfun(@(x) plot(t,x.x(1,:)), results);
% % plot(t, -max_overshoot,'--r')
% % plot(t,ss_bounds, '--k')
% % plot(t,-ss_bounds, '--k')
% % plot(t,results.exp_100.x(1,:))
% legend([dynamic_labels]); %; {'Max overhoot'; 'Steady-state error'}
% xlim([0, sim_sec]);
% ylim([-0.2, 1.2]);
% ylabel('$x$ [m]','interpreter','latex');
% xlabel('Time [s]','interpreter','latex');
% title('State Trajectories (x1)');
% subplot(2,1,2)
% hold on
% grid on
% structfun(@(x) plot(t,x.u(1,:)), results);
% legend([dynamic_labels]);
% title('State Trajectories (x2)');

% subplot(3,1,3);
% hold on
% % Extract the 'time' field from every sub-structure in 'results'
% time_evolution = structfun(@(x) x.time, results)';
% real_sim_time = ones(1,length(R_values)) * sim_sec;
% plot(R_values,time_evolution);
% plot(R_values,real_sim_time, '--r')
% xlim([min(R_values), max(R_values)]);
% legend('t', 'Real sim time');
% title('Total Runtime per N');

% subplot(3,1,3)
% hold on
% plot(t(1:M),results.exp_50.TC);
% plot(t(1:M),results.exp_50.SC);
% plot(t(1:M),results.exp_50.TC + results.exp_50.SC);
% plot(t(1:M),results.exp_100.TC);
% plot(t(1:M),results.exp_100.SC);
% plot(t(1:M),results.exp_100.TC + results.exp_100.SC);
% legend('TC50','SC50','total50','TC100','SC100','total100');
% title('Total Cost');
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


