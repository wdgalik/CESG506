%% CESG 506 HW2 - NEWTON RAPHSON w/ DISPLACEMENT CONTROL
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
u_bar = [0,0]; %Initially undeformed displacement
gamma = 0; %load factor
k_tan = EA/L1*(N1'*N1) + EA/L2*(N2'*N2); %Initial Tangent Stiffness
ev = [0,1]; % y-direction vector
k_tot = [k_tan, -Pref; -ev,0]; %Appended stiffness matrix

vert = 0; %total vertical deflection
delta = -0.01; %Vertical controlled increment of free node (m)
vert_limit = 1.2;
vert_num = vert_limit/abs(delta);

v_disp = zeros(1,vert_num);
u_disp = zeros(1,vert_num);
load_factor = zeros(1,vert_num);
Fx = zeros(1,vert_num);
Fy = zeros(1,vert_num);

tick = 0
while abs(vert) < vert_limit
    tick = tick + 1;
    load_factor(tick) = gamma;
    v_disp(tick) = delta*tick;
    R = [0;0];
    g = -delta;
    R_tilde = [R;g];  
    del_q = k_tot\R_tilde;
    delta_u = [del_q(1),del_q(2)];
    del_gamma = del_q(3);
    
    tol = 0.001;
    tock = 0;
    error = abs(delta);
    
    while error > tol
        tock = tock + 1;
        u_bar = u_bar + delta_u;
        gamma = gamma + del_gamma;
        
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
        g = ev*u_bar' - tick*delta;
        R_tilde = [R;g];
        error = norm(R_tilde);
        
        %%-----Update Tangent Stiffness-----%%
        k1 = (EA/l1_mag).*(n_1'*n_1) + (f_1/l1_mag).*(eye(2)-n_1'*n_1);
        k2 = EA/l2_mag.*(n_2'*n_2) + f_2/l2_mag.*(eye(2)-n_2'*n_2);
        k_tan = k1 + k2;
        k_tot = [k_tan, -Pref; -ev,0];
        
        %%-----Calculate New Displacements-----%%
        del_q = k_tot\R_tilde;
        delta_u = [del_q(1),del_q(2)];
        del_gamma = del_q(3);
        
        if tock > 25
            break
        else
        end
        u_disp(tick) = u_bar(1);
    end
    vert = vert + delta;
    Fx(tick) = F_int(1);
    Fy(tick) = F_int(2);
end

plot(u_disp,load_factor)
title('load factor vs. free node x-displacement')
xlabel('displacmement (m)')
ylabel('load factor, gamma')

figure
plot(v_disp,load_factor)
title('load factor vs. free node y-displacement')
xlabel('displacmement (m)')
ylabel('load factor, gamma')

figure
plot(u_disp,v_disp)
title('free node y-displacement vs. free node x-displacement')
xlabel('displacmement (m)')
ylabel('displacmement (m)')
    



