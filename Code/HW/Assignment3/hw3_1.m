%% CESG 506 HW3 - NEWTON RAPHSON ITERATION METHOD
clear; clc;

%%-----PROBLEM SPECIFIC PARAMETERS-----%%
EA = 2100; %kN
W = 5.5; %m
H = 0.5; %m
L_vec = [W,H];
L = sqrt(dot(L_vec,L_vec));
N1 = L_vec./L;
Pref = 0.3;
R = 0;
u_con = 0;
alpha = 0;
gamma_con = 0;

u_storage = [];
gamma_storage = [];

k_init = EA/L*(N1'*N1); 
k_tan = k_init(2,2);%initial stiffness

%determine delta_s & starting point for first step
dgamma = 0.25;
du = Pref*dgamma/k_tan;
delta_s = sqrt(du^2 + alpha*(Pref*dgamma)^2);
error = delta_s^2;

u_curr = u_con;
gamma_curr = gamma_con;

tick = 0;
while abs(gamma_con) < 2.0
    tick = tick + 1;
    if tick == 1
        u_curr = u_curr + du;
        gamma_curr = gamma_curr + dgamma;
    else
        plus_u = u_curr - u_con;        
        u_curr = u_curr + plus_u;
        u_con = u_con + plus_u;
        plus_gamma = gamma_curr - gamma_con;
        gamma_curr = gamma_curr + plus_gamma;
        gamma_con = gamma_con + plus_gamma;       
    end
    error = delta_s^2;
    R = 0;
    tol = 0.000001;
    tock = 0;
    while error > tol
        tock = tock + 1;
        %Find residual at updated position
        l = sqrt(W^2+(H-u_curr)^2); %m (Deformed Length)
        n = [cos(atan((H-u_curr)/W));
             sin(atan((H-u_curr)/W))];
        eps = 1/2*log(l.^2/L^2);
        f = EA.*eps;
        F_int = f.*n;
        R = gamma_curr*Pref + F_int(2);
        g = 0;
        R_tot = [R;g];
        error = norm(R_tot);
        k = (EA/l).*(n'*n) + (f/l).*(eye(2)-n'*n);
        k_tan = k(2,2);
        k_tot = [k_tan, -Pref;
                 -2*(u_curr - u_con), -2*alpha*(gamma_curr - gamma_con)];    
        q = k_tot\R_tot;
        u_curr = u_curr + q(1);
        gamma_curr = gamma_curr + q(2);
    end
    u_storage = [u_storage, u_con];
    gamma_storage = [gamma_storage, gamma_con];
end

hold on
plot(u_storage, gamma_storage*(-Pref))
title('Comparison of Load Control and ArcLength Method')
xlabel('Displacement (m)')
ylabel('Applied Force (kN)')

