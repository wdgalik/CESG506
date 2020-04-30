%% CESG 506 HW2 - NEWTON RAPHSON w/ ARC-LENGTH METHOD
clear; clc;

%%-----PROBLEM SPECIFIC PARAMETERS-----%%
EA = 2100; %kN
W1 = 5.5; %m
W2 = 4.0; %m
H = 0.5; %m
L1_vec = [W1,H];
L2_vec = [W1,H]-[W1+W2,0];
L1 = sqrt(dot(L1_vec,L1_vec));
L2 = sqrt(dot(L2_vec,L2_vec));
N1 = L1_vec./L1;
N2 = L2_vec./L2;
Pref = [0; -0.99]; %kN 

%%-----ITERATIVE NEWTON RAPHSON ALGORITHM-----%%

%%-----Initialization-----%%
u_con = [0,0]; %initial converged displacement
u_bar = u_con; %Initially undeformed displacement
u_storage = {};
u_storage = [u_storage, u_con];
gamma_con = 0;
gamma = gamma_con; %load factor
gamma_storage = [];
gamma_storage = [gamma_storage, gamma_con];

k_tan = EA/L1*(N1'*N1) + EA/L2*(N2'*N2); %Initial Tangent Stiffness

%%initialize delta_s%%
alpha = 0;
dgamma = 0.2;
du = k_tan\(Pref*dgamma);
delta_s = sqrt(dot(du,du)+alpha*dgamma^2*dot(Pref,Pref));

tick = 0;
while gamma < 2.0
    tick = tick + 1;
    %%---Converged state conditions---%%
    if tick == 1
        u_bar = u_bar + du';
        gamma = gamma + dgamma;
    else
        plus_u = u_bar - u_con;        
        u_bar = u_bar + plus_u;
        u_con = u_con + plus_u;
        u_storage = [u_storage, u_con];
        plus_gamma = gamma - gamma_con;
        gamma = gamma + plus_gamma;
        gamma_con = gamma_con + plus_gamma;  
        gamma_storage = [gamma_storage, gamma_con];
    end
    R = [0;0];
    tol = 0.0000001;
    tock = 0;
    error = delta_s^2; 
    while error > tol
        tock = tock + 1;
        %%----Find Internal Force Vector at Updated Position-----%%
        l1_vec = L1_vec + u_bar;
        l2_vec = L2_vec + u_bar;
        l1_mag = sqrt(dot(l1_vec,l1_vec));
        l2_mag = sqrt(dot(l2_vec,l2_vec));
        n_1 = l1_vec./l1_mag;
        n_2 = l2_vec./l2_mag;
        lambda_1 = l1_mag/L1;
        lambda_2 = l2_mag/L2;
        eps_1 = log(lambda_1);
        eps_2 = log(lambda_2);
        f_1 = EA*eps_1;
        f_2 = EA*eps_2;
        F_int = f_1.*n_1 + f_2.*n_2;
        
        %%-----Find Residuals-----%%
        R = -F_int' + Pref*gamma;
        g = 0; 
        R_tilde = [R;g];
        error = norm(R_tilde);
        
        %%-----Update Tangent Stiffness-----%%
        k1 = (EA/l1_mag).*(n_1'*n_1) + (f_1/l1_mag).*(eye(2)-n_1'*n_1);
        k2 = EA/l2_mag.*(n_2'*n_2) + f_2/l2_mag.*(eye(2)-n_2'*n_2);
        k_tan = k1 + k2;
        k_tot = [k_tan, -Pref; 
                -2*(u_bar-u_con),-2*alpha*(gamma-gamma_con)*dot(Pref,Pref)];
        
        %%-----Calculate New Displacements-----%%
        del_q = k_tot\R_tilde;
        u_bar = u_bar + [del_q(1),del_q(2)];
        gamma = gamma + del_q(3);
        
        if tock > 25
            break
        else
        end
    end  
end
plus_u = u_bar - u_con;        
u_con = u_con + plus_u;
u_storage = [u_storage, u_con];
plus_gamma = gamma - gamma_con;
gamma_con = gamma_con + plus_gamma;  
gamma_storage = [gamma_storage, gamma_con];

%%---Extracting Data from u_storage cells---%%
u_disp = zeros(1,tick);
v_disp = zeros(1,tick);
for i = 1:tick + 1
    u_disp(i) = u_storage{i}(1);
    v_disp(i) = u_storage{i}(2);
end

hold on
plot(v_disp, gamma_storage,'o')
xlabel('Vertical Displacement, v')
ylabel('Load Factor, gamma')
title('Load Factor vs. Vertical Displacement')
% figure

% plot(u_disp, gamma_storage,'o')
% xlabel('Horizontal Displacement, u')
% ylabel('Load Factor, gamma')
% title('Load Factor vs. HOrizontal Displacement')

    



