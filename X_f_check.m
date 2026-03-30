
Q = 1000*eye(dim_A);
R = 1*eye(dim_B);
[o,K,i] = idare(sys_d.A,sys_d.B,Q,R);

for x_1 = linspace(-10,10,100)
    for x_2 = linspace(-7,7,100)
        for x_3 = linspace(-pi/18,pi/18,100)
            for x_4 = linspace(-10,10,100)
                x_state = [x_1;x_2;x_3;x_4];
                u_k = -K*x_state;
                if(D_terminal*x_state <= c_terminal)
                   x_plus = [x_2;
                             -(81*sin(x_3)*x_4^2 + 4*u_k + (30*cos(x_3)*(u_k - (79461*sin(x_3))/100))/43)/((2430*cos(x_3)^2)/43 - 150);
                             x_4;
                             (2943*sin(x_3))/430 - (10*u_k)/1161 + (30*cos(x_3)*(81*sin(x_3)*x_4^2 + 4*u_k + (30*cos(x_3)*(u_k - (79461*sin(x_3))/100))/43))/(43*((2430*cos(x_3)^2)/43 - 150))];
                    if(D_terminal*x_plus > c_terminal)
                        x_plus
                    end
                end
            end
        end
    end
end