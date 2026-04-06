clear all;
close all;
clc;

MPC_Delft_exp_d;
MPC_Delft_exp_dn;
%%
Ts = 0.01;
sim_sec = 15;
t = 0:Ts:sim_sec;
M = sim_sec/Ts; 

results_d_raw = load('results_d.mat');
results_dn_raw = load('results_dn.mat');
results_d = results_d_raw.results;
results_dn = results_dn_raw.results;

raw_fields_d = fieldnames(results_d); 
dynamic_labels_d = strrep(raw_fields_d, 'exp_', '');
raw_fields_dn = fieldnames(results_dn); 
dynamic_labels_dn = strrep(raw_fields_dn, 'exp_', '');


figure('Position', [100, 100, 800, 530]);
subplot(2,1,1)
hold on
grid on
structfun(@(x) plot(t,x.kalman_pred(5,:)), results_d); % Plot the evolution of d_hat
xlim([0, sim_sec]);
ylim([0,0.5]);
ylabel('$\hat{d}$ [m]','interpreter','latex','FontSize', 22);
subplot(2,1,2)
hold on
grid on
structfun(@(x) plot(t,x.x_ref(1,:)), results_d);
legend('0.0', '0.1', '0.2', '0.3', '0.4','Location','best',FontSize=14); %; {'Max overhoot'; 'Steady-state error'}
ylabel('$x_{r}$ [m]','interpreter','latex','FontSize', 22)
xlabel('Time [s]','interpreter','latex','FontSize', 22);

figure('Position', [100, 100, 800, 530]);
subplot(2,1,1)
hold on
grid on
structfun(@(x) plot(t,x.kalman_pred(5,:)), results_dn); % Plot the evolution of d_hat
xlim([0, sim_sec]);
ylim([0,0.5]);
ylabel('$\hat{d}$ [m]','interpreter','latex','FontSize',22);
subplot(2,1,2)
hold on
grid on
structfun(@(x) plot(t,x.x_ref(1,:)), results_dn);
legend('0.0', '0.1', '0.2', '0.3', '0.4','Location','best',FontSize=14); %; {'Max overshoot'; 'Steady-state error'}
ylabel('$x_{r}$ [m]','interpreter','latex','FontSize', 22)
xlabel('Time [s]','interpreter','latex','FontSize', 22);
