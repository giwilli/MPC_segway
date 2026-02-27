syms m M I R x_dd phi_dd l phi phi_d u theta g x1 x2 x3 x4

eq1 = (m + M + I/R^2)*x_dd + 0.5*m*phi_dd*l*cos(phi) - 0.5*m*phi_d^2*l*sin(phi) == u/R;
eq2 = (m*(l^2/4) + theta)*phi_dd + 0.5*m*x_dd*l*cos(phi) - m*g*(l/2)*sin(phi) == -u;

phi_dd_res = solve(eq2, phi_dd);

eq1_subs = subs(eq1,phi_dd, phi_dd_res);

x_dd_res = solve(eq1_subs, x_dd);


x1_d = x2;
x2_d = subs(x_dd_res, [phi, phi_d], [x3 x4]);
x3_d = x4;
x4_d = subs(subs(phi_dd_res,x_dd_res), [phi, phi_d], [x3 x4]);

f = [x1_d; x2_d; x3_d; x4_d];

A = jacobian(f, [x1 x2 x3 x4]);

A_lin = subs(A, [x1,x2,x3,x4, u],[0,0,0,0,0])

B = jacobian(f, [u]);
B_lin = subs(B, [x1,x2,x3,x4, u],[0,0,0,0, 0])

C_lin = [1/R 0 0 0;0 0 0 1]
D_lin = [0;0];



m_n = 80 + 10; % m human + m stick
l_n = 1.8; % [m] 
R_n = 0.5/2; % [m] wheel
M_n = 40; % [kg]
I_n = 0.5*M_n*R_n^2;
theta_n = (1/3)*M_n*l_n^2;
g_n = 9.81;


T_sample = 0.01;

A_lin_s = double(subs(A_lin, [m,l,R,M,I,theta, g], [m_n, l_n,R_n, M_n, I_n, theta_n, g_n]));
B_lin_s = double(subs(B_lin, [m,l,R,M,I,theta, g], [m_n, l_n,R_n, M_n, I_n, theta_n, g_n]));
C_lin_s = double(subs(C_lin, [m,l,R,M,I,theta, g], [m_n, l_n,R_n, M_n, I_n, theta_n, g_n]));
D_lin_s = double(subs(D_lin, [m,l,R,M,I,theta, g], [m_n, l_n,R_n, M_n, I_n, theta_n, g_n]));

%%
state_space = ss(A_lin_s,B_lin_s,C_lin_s,D_lin_s)
state_space_d = c2d(state_space, T_sample, "tustin")

