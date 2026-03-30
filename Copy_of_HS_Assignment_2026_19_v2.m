% Parameters CBM problem
group = 19;   % Our group number, the task is group number dependent
xconmin = 0;  % m
xconmax = 70; % M
xcon1 = 30;   % first degradation boundary
xcon2 = 50;   % second degradation boundary
xauxmin = 0;  % 
xauxmax = 20; % N_max^gr ?
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
Ngrmax = 8;      % Again Ngrmax
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
nx = 2; % State dimension [x_con; x_aux] for each section
nd = 8; % number of deltas for each section
nz = 6; % number of auxilary z variables for each section

% x = [xcon; xaux]
Bd = zeros(nx,nd);    %just for one section (B_2 in slides)
Bd(2,7) = 1;            % d7 (equation 88 in our report)

Bz = zeros(nx,nz);      % B_3 in our report
Bz(1, 1:4) = [1 1 1 1]; % z1 + z2 + z3 + z4 (equation 84)
Bz(2, 5:6) = [1 1];      % z5 + z6  (equation 88)

% MPC data
eps = 1e-6; % machine precision, do not make it too small for numerical stability (epsilon)

% lambda = 600 + 200*rand; % generate random value between 600-800
lambda = 600;   % maintanace cost scaling
Np = 6;        % prediction horizon

% Initial state
xinit = [20+0.5*abs(group-10), 22+0.5*abs(group-10), 25+0.5*abs(group-10), 27+0.5*abs(group-10), 17+0.5*abs(group-10);
      6, 6, 7, 7, 8];

% setup the problem
x = sdpvar(nx, n, Np+1,'full');      % continuous symbolic variable of dimensions 2x5x(Np+1)
x0 = sdpvar(nx, n,'full');           % initial condition 2x5
u = intvar(nu, n, Np,'full');        % inputs 1x5xNp
delta = binvar(nd, n, Np,'full');    % binary var 8x5xNp
z = sdpvar(nz, n, Np,'full');        % auxilary var 6x5xNp

constraints = [];
constraints = [constraints, x(:,:,1) == x0]; % optimizer input parameter x0 
Jdeg = 0;
Jmaint = 0;

for k = 1:Np    
    % MLD dynamics, equality constraints allowed?
    constraints = [constraints, x(:,:,k+1) == Bd*delta(:,:,k) + Bz*z(:,:,k)];
    
    % delta1 = 1 <=> xcon >= 30  equation 47    
    constraints = [constraints, (30*delta(1,:,k) + xconmin*(1-delta(1,:,k))) <= x(1,:,k) <= ((30-eps)*(1-delta(1,:,k)) + xconmax*delta(1,:,k))];
    
    % delta2 = 1 <=> xcon >= 50  equation 47
    constraints = [constraints, (50*delta(2,:,k) + xconmin*(1-delta(2,:,k))) <= x(1,:,k) <= ((50-eps)*(1-delta(2,:,k)) + xconmax*delta(2,:,k))];
    
    % delta1delta2 = 01 cannot be true equation 11
    constraints = [constraints, (delta(2,:,k) - delta(1,:,k)) <= 0];

    % delta3 = 0 <=> xcon <= xeff equation 48
    constraints = [constraints, ((xeff+eps)*delta(3,:,k) + xconmin*(1-delta(3,:,k))) <= x(1,:,k) <= (xeff*(1-delta(3,:,k)) + xconmax*delta(3,:,k))];
        
    % delta4 = 1 <=> u >= 1 equation 53
    constraints = [constraints, (delta(4,:,k) + umin*(1-delta(4,:,k)))  <= u(1,:,k) <= ((1-eps)*(1-delta(4,:,k)) + umax*delta(4,:,k))];
    
    % delta5 = 1 <=> u >= 2 equation 54
    constraints = [constraints, (2*delta(5,:,k) + umin*(1-delta(5,:,k)))  <= u(1,:,k) <= ((2-eps)*(1-delta(5,:,k)) + umax*delta(5,:,k))];
    
    % delta4delta5 = 01 cannot be true
    constraints = [constraints, (delta(5,:,k) - delta(4,:,k)) <= 0];
    
    % delta6 = (1-d4)(1-d5) = 1 <=> u = 0 equation 77
    constraints = [constraints, (-(1-delta(4,:,k)) + delta(6,:,k)) <= 0];
    constraints = [constraints, (-(1-delta(5,:,k)) + delta(6,:,k)) <= 0];
    constraints = [constraints, ((1-delta(4,:,k)) + (1-delta(5,:,k)) - delta(6,:,k)) <= 1];
    
    % delta7 = d4(1-d5) = 1 <=> u = 1 equation 78
    constraints = [constraints, (-delta(4,:,k) + delta(7,:,k)) <= 0];
    constraints = [constraints, (-(1-delta(5,:,k)) + delta(7,:,k)) <= 0];
    constraints = [constraints, (delta(4,:,k) + (1-delta(5,:,k)) - delta(7,:,k)) <= 1];
    
    % delta8 = d4d5 = 1 <=> u = 2 equation 79
    constraints = [constraints, (-delta(4,:,k) + delta(8,:,k)) <= 0];
    constraints = [constraints, (-delta(5,:,k) + delta(8,:,k)) <= 0];
    constraints = [constraints, (delta(4,:,k) + delta(5,:,k) - delta(8,:,k)) <= 1];
    
    % z1 = d6(1-d1)(1-d2)[a1*xcon + b1] equation 80 
    constraints = [constraints, xconmin*(1-delta(1,:,k)) <= z(1,:,k) <= xconmax*(1-delta(1,:,k))];
    constraints = [constraints, xconmin*(1-delta(2,:,k)) <= z(1,:,k) <= xconmax*(1-delta(2,:,k))];
    constraints = [constraints, xconmin*delta(6,:,k) <= z(1,:,k) <= xconmax*delta(6,:,k)];
    constraints = [constraints, z(1,:,k) <= (a1*x(1,:,k) + b1 + (xconmax-xconmin)*(delta(1,:,k)+delta(2,:,k)+(1-delta(6,:,k))))];
    constraints = [constraints, z(1,:,k) >= (a1*x(1,:,k) + b1 - (xconmax-xconmin)*(delta(1,:,k)+delta(2,:,k)+(1-delta(6,:,k))))];
    
    % z2 = d6d1(1-d2)[a2*xcon + b2] equation 81 
    constraints = [constraints, xconmin*delta(1,:,k) <= z(2,:,k) <= xconmax*delta(1,:,k)];
    constraints = [constraints, xconmin*(1-delta(2,:,k)) <= z(2,:,k) <= xconmax*(1-delta(2,:,k))];
    constraints = [constraints, xconmin*delta(6,:,k) <= z(2,:,k) <= xconmax*delta(6,:,k)];
    constraints = [constraints, z(2,:,k) <= (a2*x(1,:,k) + b2 + (xconmax-xconmin)*((1-delta(1,:,k))+delta(2,:,k)+(1-delta(6,:,k))))];
    constraints = [constraints, z(2,:,k) >= (a2*x(1,:,k) + b2 - (xconmax-xconmin)*((1-delta(1,:,k))+delta(2,:,k)+(1-delta(6,:,k))))];
    
    % z3 = d6d1d2[a3*xcon + b_3] equation 82 
    constraints = [constraints, xconmin*delta(1,:,k) <= z(3,:,k) <= xconmax*delta(1,:,k)];
    constraints = [constraints, xconmin*delta(2,:,k) <= z(3,:,k) <= xconmax*delta(2,:,k)];
    constraints = [constraints, xconmin*delta(6,:,k) <= z(3,:,k) <= xconmax*delta(6,:,k)];
    constraints = [constraints, z(3,:,k) <= (a3*x(1,:,k) + b3 + (xconmax-xconmin)*((1-delta(1,:,k))+(1-delta(2,:,k))+(1-delta(6,:,k))))];
    constraints = [constraints, z(3,:,k) >= (a3*x(1,:,k) + b3 - (xconmax-xconmin)*((1-delta(1,:,k))+(1-delta(2,:,k))+(1-delta(6,:,k))))];
    
    % z4 = d7d3[psi*(xcon - xeff)] equation 83 
    constraints = [constraints, xconmin*delta(3,:,k) <= z(4,:,k) <= xconmax*delta(3,:,k)];
    constraints = [constraints, xconmin*delta(7,:,k) <= z(4,:,k) <= xconmax*delta(7,:,k)];
    constraints = [constraints, z(4,:,k) <= (psi*(x(1,:,k) - xeff) + (xconmax-xconmin)*((1-delta(3,:,k))+(1-delta(7,:,k))))];
    constraints = [constraints, z(4,:,k) >= (psi*(x(1,:,k) - xeff) - (xconmax-xconmin)*((1-delta(3,:,k))+(1-delta(7,:,k))))];

    % z5 = d6xaux equation 83
    constraints = [constraints, xauxmin*delta(6,:,k) <= z(5,:,k) <= xauxmax*delta(6,:,k)];
    constraints = [constraints, z(5,:,k) <= (x(2,:,k) + (xauxmax-xauxmin)*(delta(4,:,k)+delta(5,:,k)))];
    constraints = [constraints, z(5,:,k) >= (x(2,:,k) - (xauxmax-xauxmin)*(delta(4,:,k)+delta(5,:,k)))];
    
    % z6 = d7xaux equation 84
    constraints = [constraints, xauxmin*delta(7,:,k) <= z(6,:,k) <= xauxmax*delta(7,:,k)];
    constraints = [constraints, z(6,:,k) <= (x(2,:,k) - xauxmin*(1-delta(7,:,k)))];
    constraints = [constraints, z(6,:,k) >= (x(2,:,k) - xauxmax*(1-delta(7,:,k)))];

    % xcon(k+1) <= xmax this one comes from the need of staying below x_max
    constraints = [constraints, x(1,:,k+1) <= xmax];
    
    % xaux(k+1) <= Ngrmax this one comes from the need of staying below
    % Ngrmax
    constraints = [constraints, x(2,:,k+1) <= Ngrmax];

    Jdeg = Jdeg + sum(x(1,:,k+1)); % only sum xcon
    Jmaint = Jmaint + sum(gamma * delta([6,7,8],:,k)); % delta6,7,8 correspond to u=0,1,2 respectively
end

objective = Jdeg + lambda*Jmaint;
ops = sdpsettings('solver','glpk','verbose',1);
controller = optimizer(constraints, objective, ops, x0, {u, x});
% % controller = optimizer(constraints, objective, ops, x0, u(:,:,1));
% 
%sol = controller{xinit};
%uopt = sol{1}
%xpred = sol{2}
% % uopt = controller{xinit}

tic
diagnostics = optimize([constraints, x0 == xinit], objective, ops); % this code is only for testing, use above controller

diagnostics.problem
diagnostics.info
uopt = value(u)
xpred = value(x)
toc

%%

T = 60;
x_log = zeros(nx,n,T+1);
u_log = zeros(nu,n,T);
x_log(:,:,1) = xinit;

for t = 1:T
    disp(t)

    x0num = x_log(:,:,t);
    diagnostics = optimize([constraints, x0 == x0num], objective, ops);

    if diagnostics.problem ~= 0
        error("MPC failed at t=%d: %s", t, diagnostics.info);
    end

    uopt = value(u);
    u0 = uopt(:,:,1);
    u_log(:,:,t) = u0;

    for i = 1:n
        if u0(1,i) == 0
            if x_log(1,i,t) <= xcon1
                x_log(1,i,t+1) = a1real*x_log(1,i,t) + b1real;
            elseif x_log(1,i,t) <= xcon2
                x_log(1,i,t+1) = a2real*x_log(1,i,t) + b2real;
            else
                x_log(1,i,t+1) = a3real*x_log(1,i,t) + b3real;
            end
            x_log(2,i,t+1) = x_log(2,i,t);

        elseif u0(1,i) == 1
            if x_log(1,i,t) <= xrealeff
                x_log(1,i,t+1) = 0;
            else
                x_log(1,i,t+1) = psireal*(x_log(1,i,t) - xrealeff);
            end
            x_log(2,i,t+1) = x_log(2,i,t) + 1;

        else % u = 2
            x_log(1,i,t+1) = 0;
            x_log(2,i,t+1) = 0;
        end
    end
end

%%

subplot(3,1,1)
plot(squeeze(u_log(1,:,:))', LineStyle="none",Marker="*")
legend('Section 1','Section 2','Section 3','Section 4','Section 5')
title('Input u')


subplot(3,1,2)
plot(squeeze(x_log(1,:,:))')
legend('Section 1','Section 2','Section 3','Section 4','Section 5')
title('x_{con}')

subplot(3,1,3)
plot(squeeze(x_log(2,:,:))')
legend('Section 1','Section 2','Section 3','Section 4','Section 5')
title('x_{aux}')


