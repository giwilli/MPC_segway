% Parameters CBM problem
group = 19;   % Our group number, the task is group number dependent
xconmin = 0;  % m
xconmax = 70; % M
xcon1 = 30;   % first degradation boundary
xcon2 = 50;   % second degradation boundary
xauxmin = 0;  % 
xauxmax = 20; % 
umin = 0;     % maximum input
umax = 2;     % maximum output

% Vector A
a1 = 1.0075;  
a2 = 1.0095;
a3 = 1.01;

% Vector B
b1 = 1.65 + 0.001*abs(group-10); % group dependent
b2 = 2.07 + 0.001*abs(group-10);
b3 = 2.47 + 0.001*abs(group-10);

% These are the real system parameters, use this as our dynamics model but
% not for MPC prediction (Model-Plant mismatch)
a1real = 1.0018;
a2real = 1.016;
a3real = 1.009;
b1real = 1.6969;
b2real = 2.2175;
b3real = 2.6709;

xeff = 11.8;    % Defines squat length below which the grinding operation is
                % considered perfectly effective
xrealeff = 16;
xmax = 55;       % constraint on x^con
Ngrmax = 8;      % Ngrmax
psi = 0.96;      % Parameter representing the effectiveness of grinding
psireal = 1.008; % Real value of this parameter

xconmin_eff = xconmin - xeff;
xconmax_eff = xconmax - xeff;

% Costs for maintenance
gamma0 = 0;
gamma1 = 1;
gamma2 = 3;
gamma = [gamma0 gamma1 gamma2];

% Model data
n = 5;   % number of sections
nu = 1;  % input dimension for each section
nx = 2;  % State dimension [x_con; x_aux] for each section
nd = 8;  % number of deltas for each section
nz = 6;  % number of auxilary z variables for each section

% x = [xcon; xaux]
Bd = zeros(nx,nd);       % just for one section (B_2 in slides)
Bd(2,7) = 1;             % d7 (equation 88 in our report)

Bz = zeros(nx,nz);       
Bz(1, 1:4) = [1 1 1 1];  % z1 + z2 + z3 + z4 (equation 84)
Bz(2, 5:6) = [1 1];      % z5 + z6  (equation 88)

% MPC data
eps = 1e-6; 
lambda = 600;   
Np = 6;        

% Initial state
xinit = [20+0.5*abs(group-10), 22+0.5*abs(group-10), 25+0.5*abs(group-10), 27+0.5*abs(group-10), 17+0.5*abs(group-10);
         6, 6, 7, 7, 8];

% setup the problem
x = sdpvar(nx, n, Np+1, 'full');     
x0 = sdpvar(nx, n, 'full');          
u = intvar(nu, n, Np, 'full');       
delta = binvar(nd, n, Np, 'full');   
z = sdpvar(nz, n, Np, 'full');       

constraints = [];
constraints = [constraints, x(:,:,1) == x0];
Jdeg = 0;
Jmaint = 0;

for k = 1:Np    
    % MLD dynamics
    constraints = [constraints, x(:,:,k+1) == Bd*delta(:,:,k) + Bz*z(:,:,k)];
    
    % delta1 = 1 <=> xcon >= 30
    constraints = [constraints, ...
        (30*delta(1,:,k) + xconmin*(1-delta(1,:,k))) <= x(1,:,k) <= ...
        ((30-eps)*(1-delta(1,:,k)) + xconmax*delta(1,:,k))];
    
    % delta2 = 1 <=> xcon >= 50
    constraints = [constraints, ...
        (50*delta(2,:,k) + xconmin*(1-delta(2,:,k))) <= x(1,:,k) <= ...
        ((50-eps)*(1-delta(2,:,k)) + xconmax*delta(2,:,k))];
    
    % delta1delta2 = 01 cannot be true
    constraints = [constraints, (delta(2,:,k) - delta(1,:,k)) <= 0];

    % delta3 = 0 <=> xcon <= xeff
    constraints = [constraints, ...
        ((xeff+eps)*delta(3,:,k) + xconmin*(1-delta(3,:,k))) <= x(1,:,k) <= ...
        (xeff*(1-delta(3,:,k)) + xconmax*delta(3,:,k))];
        
    % delta4 = 1 <=> u >= 1
    constraints = [constraints, ...
        (delta(4,:,k) + umin*(1-delta(4,:,k))) <= u(1,:,k) <= ...
        ((1-eps)*(1-delta(4,:,k)) + umax*delta(4,:,k))];
    
    % delta5 = 1 <=> u >= 2
    constraints = [constraints, ...
        (2*delta(5,:,k) + umin*(1-delta(5,:,k))) <= u(1,:,k) <= ...
        ((2-eps)*(1-delta(5,:,k)) + umax*delta(5,:,k))];
    
    % delta4delta5 = 01 cannot be true
    constraints = [constraints, (delta(5,:,k) - delta(4,:,k)) <= 0];
    
    % delta6 = (1-d4)(1-d5) = 1 <=> u = 0
    constraints = [constraints, (-(1-delta(4,:,k)) + delta(6,:,k)) <= 0];
    constraints = [constraints, (-(1-delta(5,:,k)) + delta(6,:,k)) <= 0];
    constraints = [constraints, ((1-delta(4,:,k)) + (1-delta(5,:,k)) - delta(6,:,k)) <= 1];
    
    % delta7 = d4(1-d5) = 1 <=> u = 1
    constraints = [constraints, (-delta(4,:,k) + delta(7,:,k)) <= 0];
    constraints = [constraints, (-(1-delta(5,:,k)) + delta(7,:,k)) <= 0];
    constraints = [constraints, (delta(4,:,k) + (1-delta(5,:,k)) - delta(7,:,k)) <= 1];
    
    % delta8 = d4d5 = 1 <=> u = 2
    constraints = [constraints, (-delta(4,:,k) + delta(8,:,k)) <= 0];
    constraints = [constraints, (-delta(5,:,k) + delta(8,:,k)) <= 0];
    constraints = [constraints, (delta(4,:,k) + delta(5,:,k) - delta(8,:,k)) <= 1];
    
    % z1 = d6(1-d1)(1-d2)[a1*xcon + b1]
    constraints = [constraints, xconmin*(1-delta(1,:,k)) <= z(1,:,k) <= xconmax*(1-delta(1,:,k))];
    constraints = [constraints, xconmin*(1-delta(2,:,k)) <= z(1,:,k) <= xconmax*(1-delta(2,:,k))];
    constraints = [constraints, xconmin*delta(6,:,k) <= z(1,:,k) <= xconmax*delta(6,:,k)];
    constraints = [constraints, z(1,:,k) <= (a1*x(1,:,k) + b1 + (xconmax-xconmin)*(delta(1,:,k)+delta(2,:,k)+(1-delta(6,:,k))))];
    constraints = [constraints, z(1,:,k) >= (a1*x(1,:,k) + b1 - (xconmax-xconmin)*(delta(1,:,k)+delta(2,:,k)+(1-delta(6,:,k))))];
    
    % z2 = d6d1(1-d2)[a2*xcon + b2]
    constraints = [constraints, xconmin*delta(1,:,k) <= z(2,:,k) <= xconmax*delta(1,:,k)];
    constraints = [constraints, xconmin*(1-delta(2,:,k)) <= z(2,:,k) <= xconmax*(1-delta(2,:,k))];
    constraints = [constraints, xconmin*delta(6,:,k) <= z(2,:,k) <= xconmax*delta(6,:,k)];
    constraints = [constraints, z(2,:,k) <= (a2*x(1,:,k) + b2 + (xconmax-xconmin)*((1-delta(1,:,k))+delta(2,:,k)+(1-delta(6,:,k))))];
    constraints = [constraints, z(2,:,k) >= (a2*x(1,:,k) + b2 - (xconmax-xconmin)*((1-delta(1,:,k))+delta(2,:,k)+(1-delta(6,:,k))))];
    
    % z3 = d6d1d2[a3*xcon + b3]
    constraints = [constraints, xconmin*delta(1,:,k) <= z(3,:,k) <= xconmax*delta(1,:,k)];
    constraints = [constraints, xconmin*delta(2,:,k) <= z(3,:,k) <= xconmax*delta(2,:,k)];
    constraints = [constraints, xconmin*delta(6,:,k) <= z(3,:,k) <= xconmax*delta(6,:,k)];
    constraints = [constraints, z(3,:,k) <= (a3*x(1,:,k) + b3 + (xconmax-xconmin)*((1-delta(1,:,k))+(1-delta(2,:,k))+(1-delta(6,:,k))))];
    constraints = [constraints, z(3,:,k) >= (a3*x(1,:,k) + b3 - (xconmax-xconmin)*((1-delta(1,:,k))+(1-delta(2,:,k))+(1-delta(6,:,k))))];
    
    % z4 = d7d3[psi*(xcon - xeff)]
    constraints = [constraints, xconmin*delta(3,:,k) <= z(4,:,k) <= xconmax*delta(3,:,k)];
    constraints = [constraints, xconmin*delta(7,:,k) <= z(4,:,k) <= xconmax*delta(7,:,k)];
    constraints = [constraints, z(4,:,k) <= (psi*(x(1,:,k) - xeff) + (xconmax-xconmin)*((1-delta(3,:,k))+(1-delta(7,:,k))))];
    constraints = [constraints, z(4,:,k) >= (psi*(x(1,:,k) - xeff) - (xconmax-xconmin)*((1-delta(3,:,k))+(1-delta(7,:,k))))];

    % z5 = d6xaux
    constraints = [constraints, xauxmin*delta(6,:,k) <= z(5,:,k) <= xauxmax*delta(6,:,k)];
    constraints = [constraints, z(5,:,k) <= (x(2,:,k) + (xauxmax-xauxmin)*(delta(4,:,k)+delta(5,:,k)))];
    constraints = [constraints, z(5,:,k) >= (x(2,:,k) - (xauxmax-xauxmin)*(delta(4,:,k)+delta(5,:,k)))];
    
    % z6 = d7xaux
    constraints = [constraints, xauxmin*delta(7,:,k) <= z(6,:,k) <= xauxmax*delta(7,:,k)];
    constraints = [constraints, z(6,:,k) <= (x(2,:,k) - xauxmin*(1-delta(7,:,k)))];
    constraints = [constraints, z(6,:,k) >= (x(2,:,k) - xauxmax*(1-delta(7,:,k)))];

    % constraints on next state
    constraints = [constraints, x(1,:,k+1) <= xmax];
    constraints = [constraints, x(2,:,k+1) <= Ngrmax];

    Jdeg = Jdeg + sum(x(1,:,k+1));
    Jmaint = Jmaint + sum(gamma * delta([6,7,8],:,k));
end

objective = Jdeg + lambda*Jmaint;
ops = sdpsettings('solver','glpk','verbose',1);

% Precompile parametric controller
controller = optimizer(constraints, objective, ops, x0, u);

% Optional test call
sol = controller{xinit}
% uopt = sol{1}
% xpred = sol{2}

%% Closed-loop simulation using controller
T = 60;
x_log = zeros(nx,n,T+1);
u_log = zeros(nu,n,T);
x_log(:,:,1) = xinit;

for t = 1:T
    disp(t)

    x0num = x_log(:,:,t);

    % Solve MPC through optimizer object
    try
        sol = controller{x0num};
    catch ME
        error("MPC failed at t=%d: %s", t, ME.message);
    end

    % uopt = sol{1};
    % xpred = sol{2};
    uopt = sol;
    u0 = uopt(:,:,1);
    u_log(:,:,t) = u0;

    xcon_now = x_log(1,:,t);
    xaux_now = x_log(2,:,t);
    u_now    = u0(1,:); % why is this needed?
    
    xcon_next = zeros(1,n);
    xaux_next = zeros(1,n);
    
    % u = 0
    m0_1 = (u_now == 0) & (xcon_now <= xcon1);
    m0_2 = (u_now == 0) & (xcon_now > xcon1) & (xcon_now <= xcon2);
    m0_3 = (u_now == 0) & (xcon_now > xcon2);
    
    xcon_next(m0_1) = a1real*xcon_now(m0_1) + b1real;
    xcon_next(m0_2) = a2real*xcon_now(m0_2) + b2real;
    xcon_next(m0_3) = a3real*xcon_now(m0_3) + b3real;
    xaux_next(u_now == 0) = xaux_now(u_now == 0);
    
    % u = 1
    m1_1 = (u_now == 1) & (xcon_now <= xrealeff);
    m1_2 = (u_now == 1) & (xcon_now > xrealeff);
    
    xcon_next(m1_1) = 0;
    xcon_next(m1_2) = psireal*(xcon_now(m1_2) - xrealeff);
    xaux_next(u_now == 1) = xaux_now(u_now == 1) + 1;
    
    % u = 2
    m2 = (u_now == 2);
    xcon_next(m2) = 0;
    xaux_next(m2) = 0;
    
    % store
    x_log(1,:,t+1) = xcon_next;
    x_log(2,:,t+1) = xaux_next;
end

%% Plotting
ku = 0:T-1;
kx = 0:T;

figure; 
subplot(3,1,1) 
stairs(ku, squeeze(u_log(1,:,:))', LineStyle="none", Marker="*")
legend('Section 1','Section 2','Section 3','Section 4','Section 5', Location='bestoutside')
title('Input u')
xlabel('Time step [k]')
ylabel('maintenance intervention u')
xlim([0 T])
grid on;

subplot(3,1,2)
stairs(kx, squeeze(x_log(1,:,:))')
legend('Section 1','Section 2','Section 3','Section 4','Section 5', Location='bestoutside')
title('x_{con}')
xlabel('Time step [k]')
ylabel('condition x_{con} [mm]')
xlim([0 T])
grid on;

subplot(3,1,3)
stairs(kx, squeeze(x_log(2,:,:))')
legend('Section 1','Section 2','Section 3','Section 4','Section 5', Location='bestoutside')
title('x_{aux}')
xlabel('Time step [k]')
ylabel('number of grinding x_{aux}')
xlim([0 T])
grid on;