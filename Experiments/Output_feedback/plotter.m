clear all;
close all;
clc;

MPC_Delft_fullstate_exp_Q;
MPC_Delft_fullstate_exp_R;
%%
results_Q_raw = load('results_Q.mat');
results_R_raw = load('results_R.mat');
results_Q = results_Q_raw.results;
results_R = results_R_raw.results;
sim_sec = 10;

raw_fields_Q = fieldnames(results_Q); 
dynamic_labels_Q = strrep(raw_fields_Q, 'exp_', '');
raw_fields_R = fieldnames(results_R); 
dynamic_labels_R = strrep(raw_fields_R, 'exp_', '');

ss_bounds = ones(1,M+1) * 0.02;
max_overshoot = ones(1,M+1)* 0.1;

figure('Position', [100, 100, 800, 530]);
subplot(2,1,1)
hold on
grid on
structfun(@(x) plot(t,x.x(1,:)), results_Q);
legend([dynamic_labels_Q],FontSize=14); %; {'Max overhoot'; 'Steady-state error'}
xlim([0, sim_sec]);
ylim([-0.2, 1.2]);
ylabel('$x$ [m]','interpreter','latex','FontSize', 22);
subplot(2,1,2)
hold on
grid on
structfun(@(x) plot(t,x.x(1,:)), results_R);
legend('0.000001','0.00001','0,0001','0.001', '0.01',FontSize=14); %; {'Max overhoot'; 'Steady-state error'}
xlim([0, sim_sec]);
ylim([-0.2, 1.2]);
ylabel('$x$ [m]','interpreter','latex','FontSize', 22)
xlabel('Time [s]','interpreter','latex','FontSize', 22);

figure('Position', [100, 100, 800, 265]);
subplot(1,1,1)
hold on
grid on
structfun(@(x) plot(t,x.u(1,:)), results_Q);
legend([dynamic_labels_Q],FontSize=14); %; {'Max overhoot'; 'Steady-state error'}
xlim([0, sim_sec]);
ylabel('$u$ [Nm]','interpreter','latex','FontSize', 22);
xlabel('Time [s]','interpreter','latex','FontSize', 22);

% subplot(2,1,2)
% hold on
% grid on
% structfun(@(x) plot(t,x.u(1,:)), results_R);
% legend('0.000001','0.00001','0,0001','0.001', '0.01'); %; {'Max overhoot'; 'Steady-state error'}
% xlim([0, sim_sec]);
% ylabel('$x$ [m]','interpreter','latex')
% xlabel('Time [s]','interpreter','latex');