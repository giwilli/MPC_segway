A = A_lin_s;
B = B_lin_s;
C = C_lin_s;
D = D_lin_s;

Ts = 0.01;
sys_d = c2d(ss(A,B, C, D),Ts,'tustin');

rank(ctrb(sys_d.A, sys_d.B))
%%

N = 100;
dim_A = size(A,1);
dim_B = size(B,2);
dim_C = size(C,1);
Q = 1000*eye(dim_A);
R = 1*eye(dim_B);

[P,K,L] = idare(A,B,Q,R);


P = 1000*P;
x0 = [-0.1;-0.1;-0.1;-0.1];

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
c = [Inf; 7; pi/18; Inf];
u_bound = 42;

D_tilde = [D*sys_d.A;-D*sys_d.A; zeros(1,4); zeros(1,4)];
E_tilde = [D*sys_d.B;-D*sys_d.B; 1; -1];
b_tilde =[c;c;u_bound;u_bound];
D_bar = kron(eye(N),D_tilde);
E_bar = kron(eye(N),E_tilde);
b_bar = repmat(b_tilde,[N,1]);
T_tilde = T(1:end-dim_A,:);
S_tilde = S(1:end-dim_A,:);

G = D_bar*S_tilde + E_bar;
g = b_bar - D_bar*T_tilde*x0;

sys =  LTISystem('A', sys_d.A, 'B', sys_d.B, 'C', sys_d.C, 'D', sys_d.D);
sys.x.min = -c;
sys.x.max = c;
sys.u.min = -u_bound;
sys.u.max = u_bound;

Xf = sys.invariantSet();
Xf.plot()




%%
%u = quadprog(H,h, G,g);


M = 2000;
x = zeros(dim_A, (M+1));
x(:,1) = x0;

for i = 1:M
    x0 = x(:,i);
    h = S.'*Q_bar*T*x0;
    g = b_bar - D_bar*T_tilde*x0;
    u = quadprog(H,h, G,g);
    x(:,i+1) = sys_d.A*x(:,i) + sys_d.B*u(1);
end
%%



%%

eig(sys_d.A)
