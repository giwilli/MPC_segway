A = A_lin_s;
B = B_lin_s;
C = C_lin_s;
D = D_lin_s;




Ts = 0.01;
sys_d = c2d(ss(A,B, C, D),Ts,'tustin');

rank(ctrb(sys_d.A, sys_d.B))
%%

global mptOptions
mptOptions.verbose = 1;

%%

N = 70;
dim_A = size(A,1);
dim_B = size(B,2);
dim_C = size(C,1);
Q = 1000*eye(dim_A);
R = 1*eye(dim_B);
x0 = [-0.1;-0.1;-0.1;-0.1];

%%
c = [Inf; 7; pi/18; Inf];
u_bound = 42;

model = LTISystem(sys_d);
model.x.min = -c;
model.x.max = c;
model.u.min = -u_bound;
model.u.max = u_bound;

model.x.penalty = QuadFunction(Q);
model.u.penalty = QuadFunction(R);



Tset = model.LQRSet;
P = model.LQRPenalty.weight;


T = zeros(dim_A*(N+1),dim_A);
S = zeros(dim_A*(N+1), dim_B*N);


for i = 1:N+1
    T((i-1)*dim_A + 1 : i*dim_A, 1:dim_A) = sys_d.A^(i-1);
    if i <= N
        S(:,end-i*dim_B + 1:end -(i-1)*dim_B) = circshift(T,-dim_A*i,1)*sys_d.B;
    end
end

Q_bar = blkdiag(kron(eye(N),Q), P);
R_bar = kron(eye(N),R);

H = S.'*Q_bar*S + R_bar;
H = (H+H')/2;

h = S.'*Q_bar*T*x0;

D = diag([0 1 1 0]);

D_tilde = [D*sys_d.A;-D*sys_d.A; zeros(1,dim_A); zeros(1,dim_A)];
E_tilde = [D*sys_d.B;-D*sys_d.B; 1; -1];
b_tilde =[c;c;u_bound;u_bound];
D_bar_temp = kron(eye(N),D_tilde);
E_bar_temp = kron(eye(N),E_tilde);
b_bar_temp = repmat(b_tilde,[N,1]);
T_tilde = T(1:end-dim_A,:);
S_tilde = S(1:end-dim_A,:);

%%


D_terminal = Tset.A;
c_terminal = Tset.b;

D_tilde_term = [D_terminal*sys_d.A;-D_terminal*sys_d.A; zeros(1,dim_A); zeros(1,dim_A)];
E_tilde_term = [D_terminal*sys_d.B;-D_terminal*sys_d.B; 0; 0];
b_tilde_term = [c_terminal;c_terminal;0;0];

tmp = zeros(1,N);
tmp(:,end) = 1;
D_bar_term_temp = kron(tmp,D_tilde_term);
E_bar_term_temp = kron(tmp,E_tilde_term);
b_bar_term_temp = b_tilde_term;

D_bar = [D_bar_temp;D_bar_term_temp];
E_bar = [E_bar_temp;E_bar_term_temp];
b_bar = [b_bar_temp;b_bar_term_temp];

%%

G = D_bar*S_tilde + E_bar;
g = b_bar - D_bar*T_tilde*x0;

%%
%u = quadprog(H,h, G,g);


M = 1000;
x = zeros(dim_A, (M+1));
x(:,1) = x0;
u_log = zeros(dim_B, M+1);
u_log(:,1) = 0;

for i = 1:M
    disp(i);
    x0 = x(:,i);
    h = S.'*Q_bar*T*x0;
    g = b_bar - D_bar*T_tilde*x0;
    u = quadprog(H,h, G,g);
    u_log(:,i) = u;
    x(:,i+1) = sys_d.A*x(:,i) + sys_d.B*u(1);
end
%% LQR for reference
[K,S,P] = lqr(sys_d,Q,R);
M = 1000;
x_lqr = zeros(dim_A, (M+1));
x_lqr(:,1) = x0;
u_log_lqr = zeros(dim_B, M+1);
u_log_lqr(:,1) = 0;

for i = 1:M
    %disp(i);
    u = -K*x(:,i);
    x_lqr(:,i+1) = (sys_d.A*x_lqr(:,i) + sys_d.B*u);
    u_log_lqr(:,i) = u;
end
%%
hold on
t = 0:0.01:M*0.01;
plot(t,x)
plot(t, u_log)

plot(t,x_lqr)
plot(t,u_log_lqr)
