%% Cleanup and initialize
clear all
clc;
tbxmanager restorepath;
mpt_init;

%% Global Solver option
global MPTOPTIONS
MPTOPTIONS.modules.geometry.sets.Polyhedron.projection.method = 'mplp';
%% Create the discretized system
clear all;
close all;
clc;
disp('RESET');
Equations;
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

    % Running for different Q
Q_values = [1,10,100,1000];
N = 45;

    % Initialize a structure to hold your results
results = struct();

for j = 1:length(Q_values)
    Q_val = Q_values(j);
    disp('Current Q is: ');
    disp(Q_val);
    
    Q = Q_val*eye(dim_A);
    R = 1*eye(dim_B);
    
    c = [100; 7; pi/18; 100];
    u_bound = 42;
    
        % Constraints definition and Terminal Set
    
    model = LTISystem(sys_d);
    model.x.min = -c;
    model.x.max = c;
    model.u.min = -u_bound;
    model.u.max = u_bound;
    
    model.x.penalty = QuadFunction(Q);
    model.u.penalty = QuadFunction(R);
    P = model.LQRPenalty.weight;
    load_TSet = false;
    if load_TSet
        Tset_Aload = load("data/Tset_A_Q1000R1.mat");
        Tset_bload = load("data/Tset_b_Q1000R1.mat");
    
        Tset_A = Tset_Aload.Tset_A;
        Tset_b = Tset_bload.Tset_b;
    else
        Tset = model.LQRSet;
        Tset_A = Tset.A;
        Tset_b = Tset.b;
    end
    D_terminal = Tset_A;
    c_terminal = Tset_b;
        % Compact MPC Formulation
    
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
    
        % Objective function weights (Compact Form) 
    %P = load("P.mat");
    %P = P.P;
    
    Q_bar = blkdiag(kron(eye(N),Q), P);
    R_bar = kron(eye(N),R);
    
    H = S.'*Q_bar*S + R_bar;
    H = (H+H')/2;
    
    %h = S.'*Q_bar*T*x0;

        % Terminal Constraint Formulation
    
    % 
    
    D_tilde_term = [D_terminal*sys_d.A];   %;-D_terminal*sys_d.A; zeros(1,dim_A); zeros(1,dim_A)];
    E_tilde_term = [D_terminal*sys_d.B];   %;-D_terminal*sys_d.B; 0; 0];
    b_tilde_term = [c_terminal];           %;c_terminal;0;0];
    
    tmp = zeros(1,N);
    tmp(:,end) = 1;
    D_bar_term_temp = kron(tmp,D_tilde_term);
    E_bar_term_temp = kron(tmp,E_tilde_term);
    b_bar_term_temp = b_tilde_term;
    
        % Constraint Concatination
    x0 = [0;0;0;0];
    
    D_bar = [D_bar_temp;D_bar_term_temp];
    E_bar = [E_bar_temp;E_bar_term_temp];
    b_bar = [b_bar_temp;b_bar_term_temp];
    
    G = D_bar*S_tilde + E_bar;
    
    g = b_bar - D_bar*T_tilde*x0;
    
        % Closed Loop global paramters
    sim_sec = 30;
    t = 0:Ts:sim_sec;
    M = sim_sec/Ts;
    y_ref_final = 0;
    y_reg = zeros(1,M+1);
    y_constant = ones(1,M+1)* y_ref_final;
    y_square = square(pi*t/10);
    y_sine = sin(t);
    y_linear = linspace(0,y_ref_final,M);
    y_step = (t>4);
    
        %Reference given to the controller
    y_ref = y_constant;%[linspace(0,y_ref_final,M)];

        % Closed Loop MPC
    x_mpc = zeros(dim_A, (M+1));
    y_mpc = zeros(dim_C, (M+1));
    x_mpc(:,1) = x0;
    y_mpc(:,1) = 0;
    u_mpc_log = zeros(dim_B, M+1);
    u_mpc_log(:,1) = 0;
    n_d = 1;
    
    kalman_log = zeros(dim_A+n_d, (M+1));
    x_ref_normlog = zeros(dim_A, (M+1));
    
    d = 0.0;
    
    q_std = 0.00001;
    pos_std = 0.002;
    vel_std = 0.002;
    angv_std = 0.02;
    
    P_Kalm =  (1e-3)*eye(dim_A+n_d);
    Q_Kalm = q_std^2*eye(dim_A+n_d);
    R_Kalm = diag([pos_std^2 vel_std^2 angv_std^2]);
    
    w = sqrt(Q_Kalm)*randn(dim_A+n_d,M+1)* 0;
    v = sqrt(R_Kalm)*randn(dim_C,M+1)*0;
    
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
        c_terminal = Tset_b + D_terminal*x_ref;
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
    % Create the dynamic field name (e.g., 'x_10')
    expName = sprintf('exp_%d', Q_val); 
    
    % Save x and u inside that specific experiment's field
    results.(expName).x = x_mpc;
    results.(expName).u = u_mpc_log;
end
%% Plot the different horizons
raw_fields = fieldnames(results); 
dynamic_labels = strrep(raw_fields, 'exp_', '');

subplot(2,1,1)
hold on
grid on
structfun(@(x) plot(t,x.x(1,:)), results);
plot(t,y_ref, '--k')
% plot(t,results.exp_100.x(1,:))
legend([dynamic_labels; {'Reference'}]);
title('State Trajectories (x1)');
subplot(2,1,2)
hold on
grid on
structfun(@(x) plot(t,x.u(1,:)), results);
legend([dynamic_labels]);
title('Input (u)');