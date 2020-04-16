%% CESG 506 HW1 - SNAP-THROUGH LOAD
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

%%-----DISPLACEMENT CONTROL (SOLVING FREE NODE TRAJECTORY-----%%
v = linspace(0,2.5*H,100); %Controlling downwards displacement
u_storage = zeros(1,100); %Will be used to store indices of u-values that satisfy equilibrium. One associated with each v disp. 
F_int_x = zeros(100,2000);
F_int_x(:,1) = 1; %Assign any positive value to the first column. Used to give initial sign value of F_x in horizontal disp loops
for i = 1:100 %loop through y displacements
    for j = 2:2000  %loop through x-displacements
        u = j/120000; 
        u_bar = [-u,-v(i)]; %displacement of free node. Down due to applied force. Left to maintain x-direction equil'b
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
        F_int_x(i,j) = F_int(:,1); 
        
        if F_int_x(i,j)*F_int_x(i,j-1) < 0
            u_storage(i) = -u; %Finding and storing the u displacements associated with v(i) 
        else
        end
        
    end
end

%%-----SOLVING AND PLOTTING INTERNAL FORCE (y) FOR FREE NODE TRAJECTORY-----%%
FF_int_y = zeros(1,100);
for ii = 1:100
    vv = v(ii); %has opposite sign
    uu = u_storage(ii); %already has correct sign
    uu_bar = [uu, -vv]; %signs both correct
    el1_vec = L1_vec + uu_bar;
    el2_vec = L2_vec + uu_bar;
    el1_mag = sqrt(dot(el1_vec,el1_vec));
    el2_mag = sqrt(dot(el2_vec,el2_vec));
    nn_1 = el1_vec./el1_mag;
    nn_2 = el2_vec./el2_mag;
    llambda_1 = el1_mag/L1;
    llambda_2 = el2_mag/L2;
    eeps_1 = log(llambda_1);
    eeps_2 = log(llambda_2);
    ff_1 = EA*eeps_1;
    ff_2 = EA*eeps_2;
    FF_int = ff_1.*nn_1 + ff_2.*nn_2;
    FF_int_y(ii) = FF_int(:,2);
end

plot(v,FF_int_y)
xlabel('Vertical Disp, v (m)')
ylabel('Applied Force,P (kN)')
title('Non-Linear Load-Displacement Curve')

        
        