%% Cleanup and initialize
clc;
tbxmanager restorepath;
mpt_init;

%% Global Solver option
global MPTOPTIONS
MPTOPTIONS.modules.geometry.sets.Polyhedron.projection.method = 'mplp';
global mptOptions
mptOptions.verbose = 1;
mptOptions.infbound = 500; % Increase from the default 100

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

rank(ctrb(sys_d.A, sys_d.B))
rank(obsv(sys_d.A, sys_d.C))

%% Problem Fundamentals

N = 50;
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

% P = load("P.mat");
% P = P.P;


Q_bar = blkdiag(kron(eye(N),Q), P);
R_bar = kron(eye(N),R);

H = S.'*Q_bar*S + R_bar;
H = (H+H')/2;

%h = S.'*Q_bar*T*x0;
%% Terminal Constraint Formulation

Tset = model.LQRSet;
Tset_A = Tset.A;
Tset_b = Tset.b;

D_terminal = Tset.A;
c_terminal = Tset.b;

% Tset_A = load("TS_A.mat");
% Tset_b = load("TS_b.mat");
% 
% D_terminal = Tset_A.TS_A;
% c_terminal = Tset_b.TS_b;

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

% x_mpc = zeros(dim_A, (M+1));
% x_mpc(:,1) = x0;
% u_mpc_log = zeros(dim_B, M+1);
% u_mpc_log(:,1) = 0;
% 
% for i = 1:M
%     disp(i);
%     x0 = x_mpc(:,i);
%     h = S.'*Q_bar*T*x0;
%     g = b_bar - D_bar*T_tilde*x0;
%     u = quadprog(H,h, G,g);
%     u_mpc_log(:,i) = u(1,:);
%     x_mpc(:,i+1) = sys_d.A*x_mpc(:,i) + sys_d.B*u(1);
% end
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

%%
bnd_r = 6.5;
bnd_s = pi/18;
r = 1;
res = 25;
mat_plot = zeros(res,res);
disp('Calculating: ')
% loop through different initial conditions 
for initial_r = linspace(-bnd_r,bnd_r,res)
    disp(r)
    s = 1;
    for initial_s = linspace(-bnd_s,bnd_s,res)
        x0 = [initial_r 0 initial_s 0]';
        h = S.'*Q_bar*T*x0;
        g = b_bar - D_bar*T_tilde*x0;
        options = optimoptions('quadprog', 'Display', 'none'); % Keep command window clean
        [u, fval, exitflag, output, lambda] = quadprog(H, h, G, g, [], [], [], [], [], options);
        % Listen to quadprogs warning
        if exitflag == 1
            mat_plot(r,s) = 1;
        end
        s = s+1;
    end
    r = r+1;
end

fprintf('\tdone!\n');
%%
% plot the initial conditions which were steered towards the terminal set
% X_f within N steps

fprintf('\tplotting results ... \n');
figure(1);
hold on;

r = 1;
m_size = 11;
for r_plot = linspace(-bnd_r,bnd_r,res)
    s = 1;
    disp(r)
    for s_plot = linspace(-bnd_s,bnd_s,res)
        if (mat_plot(r,s) == 1)
            plot(r_plot,s_plot,'sg','MarkerFaceColor','g','MarkerEdgeColor','g','MarkerSize',m_size);
        else
            plot(r_plot,s_plot,'sr','MarkerFaceColor','w','MarkerEdgeColor','w','MarkerSize',m_size);
        end
        s = s+1;
    end
    r = r+1;
end
xlabel('x [m]', 'interpreter','latex');
ylabel('$\varphi$ [rad]','interpreter','latex');
grid on;
xlim([-1.05*bnd_r 1.05*bnd_r]);
ylim([-1.05*bnd_s 1.05*bnd_s]);
fprintf('\tdone! \n');

%%
% Plot the convex hull of the initial conditions that can be steered
% to the terminal set X_f within N steps

fprintf('\tplotting results ... \n');
figure(2); clf;
hold on;

% Grid corresponding to mat_plot
r_vals = linspace(-bnd_r, bnd_r, res);
s_vals = linspace(-bnd_s, bnd_s, res);

% Build coordinate matrices so indices match mat_plot(r,s)
[R, S] = ndgrid(r_vals, s_vals);

% Extract all positive points
idx_pos = (mat_plot == 1);
x_pos = R(idx_pos);
y_pos = S(idx_pos);

if ~isempty(x_pos)
    if numel(x_pos) >= 3
        % Compute convex hull
        k = convhull(x_pos, y_pos);

        % Filled convex hull
        fill(x_pos(k), y_pos(k), 'g', ...
            'FaceAlpha', 0.3, ...
            'EdgeColor', 'g', ...
            'LineWidth', 1.5);

        % Optional: plot hull vertices / boundary
        plot(x_pos(k), y_pos(k), 'g-', 'LineWidth', 1.5);
    else
        % Fallback if too few points for a hull
        plot(x_pos, y_pos, 'sg', ...
            'MarkerFaceColor', 'g', ...
            'MarkerEdgeColor', 'g', ...
            'MarkerSize', 8);
    end
else
    warning('No positive points found in mat_plot.');
end

xlabel('x [m]', 'Interpreter', 'latex');
ylabel('$\varphi$ [rad]', 'Interpreter', 'latex');
grid on;
%axis equal;
xlim([-1.05*bnd_r, 1.05*bnd_r]);
ylim([-1.05*bnd_s, 1.05*bnd_s]);

fprintf('\tdone! \n');
%%
% Plot convex hull of feasible points, and fill outside region in red

fprintf('\tplotting results ... \n');
figure(2); clf;
hold on;

% Grid corresponding to mat_plot
r_vals = linspace(-bnd_r, bnd_r, res);
s_vals = linspace(-bnd_s, bnd_s, res);
[R, S] = ndgrid(r_vals, s_vals);

% Extract feasible points
idx_pos = (mat_plot == 1);
x_pos = R(idx_pos);
y_pos = S(idx_pos);

% Plot limits
xL = [-1.05*bnd_r, 1.05*bnd_r];
yL = [-1.05*bnd_s, 1.05*bnd_s];

% Fill whole domain red first
fill([xL(1) xL(2) xL(2) xL(1)], ...
     [yL(1) yL(1) yL(2) yL(2)], ...
     'r', 'EdgeColor', 'none','FaceAlpha', 0.40);

if ~isempty(x_pos) && numel(x_pos) >= 3
    % Convex hull of feasible set
    k = convhull(x_pos, y_pos);

    % Fill feasible region on top
    fill(x_pos(k), y_pos(k), 'g', ...
        'FaceAlpha', 0.50, ...
        'EdgeColor', 'g', ...
        'LineWidth', 1.5);

    plot(x_pos(k), y_pos(k), 'k-', 'LineWidth', 1.5, 'Color', 'g');
elseif ~isempty(x_pos)
    plot(x_pos, y_pos, 'sg', ...
        'MarkerFaceColor', 'g', ...
        'MarkerEdgeColor', 'g', ...
        'MarkerSize', 8);
else
    warning('No positive points found in mat_plot.');
end

xlabel('x [m]', 'Interpreter', 'latex');
ylabel('$\varphi$ [rad]', 'Interpreter', 'latex');
grid on;
xlim(xL);
ylim(yL);

fprintf('\tdone! \n');

%% Xf
bnd_r = 6.5;
bnd_s = pi/18;
r = 1;
res = 25;
mat_plot_f = zeros(res,res);
disp('Calculating: ')
% loop through different initial conditions 
for initial_r = linspace(-bnd_r,bnd_r,res)
    disp(r)
    s = 1;
    for initial_s = linspace(-bnd_s,bnd_s,res)
        x0 = [initial_r 0 initial_s 0]';
        if Tset_A * x0 <= Tset_b
            mat_plot_f(r,s) = 1;
        end
        s = s+1;
    end
    r = r+1;
end

fprintf('\tdone!\n');

%% Plot Xf
% plot the initial conditions which were steered towards the terminal set
% X_f within N steps

fprintf('\tplotting results ... \n');
figure(2);
hold on;

r = 1;
m_size = 11;
for r_plot = linspace(-bnd_r,bnd_r,res)
    s = 1;
    disp(r)
    for s_plot = linspace(-bnd_s,bnd_s,res)
        if (mat_plot_f(r,s) == 1)
            plot(r_plot,s_plot,'sg','MarkerFaceColor','g','MarkerEdgeColor','g','MarkerSize',m_size);
        else
            plot(r_plot,s_plot,'sr','MarkerFaceColor','w','MarkerEdgeColor','w','MarkerSize',m_size);
        end
        s = s+1;
    end
    r = r+1;
end
xlabel('x [m]', 'interpreter','latex');
ylabel('$\varphi$ [rad]','interpreter','latex');
grid on;
xlim([-1.05*bnd_r 1.05*bnd_r]);
ylim([-1.05*bnd_s 1.05*bnd_s]);
fprintf('\tdone! \n');

