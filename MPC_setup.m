clc;
clear;

%% Initialize system dynamics

[Kf, m_h, m_w, m_f, m_b, Lh, La, Lw, g, K_EC_T, K_EC_P, K_EC_E ] = setup_heli_3d_configuration();
HELI3D_ABCD_eqns;

Ts = 0.05;
sys_d = c2d(ss(A,B, C, D),Ts,'tustin');

%% Initilize horrizon length
N = 15;

% Pick the cost matrices
Q = [1 0 0; 0 1 0;0 0 1];
R = [0.01 0;0 0.01];
S = 100*eye(size(C,1));

% Set constraints on states and inputs

x_max = [pi/3; pi/3; pi; pi/4; pi/4; pi/10];
x_min = [-pi/3; -pi/3; -pi; -pi/4; -pi/4; -pi/10];

u_max = [24; 24];
u_min = [-24; -24];

% Augment the system by introducing input increments

A_tilde = [sys_d.A sys_d.B; zeros(size(B,2), size(A,2)) eye(size(B,2))];
B_tilde = [sys_d.B;eye(size(B,2))];
C_tilde = [sys_d.C zeros(size(C,1), size(B,2))];
D_tilde = zeros(size(C,1), size(B,2));

% Define building blocks of MPC

state_vec_size = size(A_tilde,1);
input_size = size(B_tilde,2);


blocks_Q_tilde = repmat({Q}, 1, N-1);   % 1×N cell array, each entry is A
blocks_C_tilde = repmat({C_tilde},1,N);
blocks_R_tilde = repmat({R}, 1, N);
blocks_M_tilde = repmat({Q*C_tilde},1,N-1);
blocks_K_input = repmat({[zeros(input_size,state_vec_size-2) eye(input_size)]}, 1, N);
blocks_K_states = repmat({[eye(state_vec_size-2) zeros(state_vec_size-2, input_size)]}, 1, N+1);
blocks_I_tilde = repmat({eye(input_size)},1, N);

Q_d_tilde = blkdiag(blocks_Q_tilde{:}, S);
C_d_tilde = blkdiag(blocks_C_tilde{:});
R_d_tilde = blkdiag(blocks_R_tilde{:});
M_d_tilde = blkdiag(blocks_M_tilde{:}, S*C_tilde);
K_input = blkdiag(blocks_K_input{:});
K_states = blkdiag(blocks_K_states{:});
I_d_tilde = blkdiag(blocks_I_tilde{:});

Kapa = zeros(N*state_vec_size,N*input_size);
Psi = zeros(N*state_vec_size,state_vec_size);

i_index = 1;
for i = 1:N
    j_index = 1;
    for j = i-1:-1:0
        Kapa(i_index:i_index + state_vec_size - 1, j_index:j_index+input_size-1) = (A_tilde^j)*B_tilde;
        j_index = j_index + input_size;
    end
    Psi(i_index:i_index+state_vec_size-1, :) = A_tilde^i;
    i_index = i_index + state_vec_size;
end

Psi_extended_input = [eye(state_vec_size); Psi(1:(N-1)*state_vec_size,:)];
L_input = K_input*Psi_extended_input;

Psi_extended_states = [eye(state_vec_size); Psi];
L_states = K_states*Psi_extended_states;

G_input = I_d_tilde + K_input*[zeros(state_vec_size, input_size*N);Kapa(1:state_vec_size*(N-1),:)];
G_states = K_states*[zeros(state_vec_size, N*input_size);Kapa];

H = (Kapa.'*C_d_tilde.'*Q_d_tilde*C_d_tilde*Kapa + R_d_tilde);
F = [Psi.'*C_d_tilde.'*Q_d_tilde.'*C_d_tilde*Kapa;-M_d_tilde*Kapa];
G = [G_input;G_states;-G_input;-G_states];

u_min_hat = repmat(u_min, N, 1);
u_max_hat = repmat(u_max, N, 1); 
x_min_hat = repmat(x_min, N+1, 1);
x_max_hat = repmat(x_max, N+1, 1);

% Create a bus object for passing matrices into a matlab function block

elems(1) = Simulink.BusElement;
elems(1).Name = 'H';
elems(1).Dimensions = size(H);

elems(2) = Simulink.BusElement;
elems(2).Name = 'F';
elems(2).Dimensions = size(F);

elems(3) = Simulink.BusElement;
elems(3).Name = 'G';
elems(3).Dimensions = size(G);

elems(4) = Simulink.BusElement;
elems(4).Name = 'L_input';
elems(4).Dimensions = size(L_input);

elems(5) = Simulink.BusElement;
elems(5).Name = 'L_states';
elems(5).Dimensions = size(L_states);

elems(6) = Simulink.BusElement;
elems(6).Name = 'u_min_hat';
elems(6).Dimensions = size(u_min_hat);

elems(7) = Simulink.BusElement;
elems(7).Name = 'u_max_hat';
elems(7).Dimensions = size(u_max_hat);

elems(8) = Simulink.BusElement;
elems(8).Name = 'x_min_hat';
elems(8).Dimensions = size(x_min_hat);

elems(9) = Simulink.BusElement;
elems(9).Name = 'x_max_hat';
elems(9).Dimensions = size(x_max_hat);


myBus = Simulink.Bus;
myBus.Elements = elems;

%% The Optimization task tryout

x_t = [1;1;1;1;1;1;1;1];
y_des = repmat([2;2;2],N,1);

P = [u_max_hat - L_input * x_t;
     x_max_hat - L_states * x_t;
     L_input*x_t - u_min_hat;
     L_states*x_t - x_min_hat];

f = ([x_t.' y_des.']*F).';


H = 0.5*(H+H');
delta_u_opt = quadprog(H,f,G,P);



