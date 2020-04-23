%% CESG 506 HW2 - DISPLACEMENT CONTROL FRAME SYSTEM
clear; clc;
%%-----PROBLEM SPECIFIC PARAMETERS-----%%
EA = 2100; %kN
Pref = [0; 0; -0.99; 0; 0; 0]; %kN 
ek = [0,0,1,0,0,0];

%%-----DEFINE SYSTEM LAYOUT-----%%
layout = [0, 0, 0, 0; %defines original nodal layout (x,y,z,fixity)
          9.5, 0, 0, 0; %row m is node_m coordinates
          0, -5, 0, 0;
          9.5, -5, 0, 0;
          5.5, -1.25, 0.5, 1;
          5.5, -3.75, 0.5, 1];

logic = [1,1,5; %defines element connectivity (element#,nodei,nodej)
         2,1,6;
         3,2,5;
         4,3,6;
         5,4,5;
         6,4,6;
         7,5,6];
     
for i = 1:length(logic)
    Length{i} = layout(logic(i,3),(1:3))-layout(logic(i,2),(1:3));
end

%%-----INITIALIZATION-----%%
disp5 = [0,0,0];
disp6 = [0,0,0];
u_bar = [disp5,disp6];
diff = disp6 - disp5;

%%-----Get Initial Element Stiffnesses-----%%
for m = 1:length(logic) %m = element number
    if logic(m,3)==5 
        [F{m},k{m}] = stiffness(EA, Length{m}, disp5, 3);
    elseif logic(m,2)~=5
        [F{m},k{m}] = stiffness(EA, Length{m}, disp6, 3);
    else
        [F{m},k{m}] = stiffness(EA, Length{m}, diff, 3);
    end 
end

kff1 = k{1} + k{5} + k{3} + k{7}; 
kff2 = -k{7}; %Should use fixity values from layout instead of hardcoding
kff3 = -k{7};
kff4 = k{2} + k{4} + k{6} + k{7};
k_tan = [kff1,kff2;kff3,kff4];

k_tot = [k_tan,-Pref;-ek,0]; %Appended stiffness matrix

%%----ITERATIVE NEWTON-RAPHSON ALGORITHM-----%%
vert = 0; %total vertical deflection (z-dir)(controlled variable)
delta = -0.005; %Vertical controlled increment of node 5 (m)

u_disp5 = [];
v_disp5 = [];
w_disp5 = [];
u_disp6 = [];
v_disp6 = [];
w_disp6 = [];

load_factor = [];
gamma = 0;
F36 = 0;
F46 = 0;

tick = 0;
while abs(vert) < 1.15
    tick = tick + 1;
    load_factor(tick) = gamma;
    w_disp5(tick) = delta*tick;
    R = [0;0;0;0;0;0];
    g = -delta;
    R_tilde = [R;g];  
    del_q = k_tot\R_tilde;
    delta_u = [del_q(1),del_q(2),del_q(3),del_q(4),del_q(5),del_q(6)];
    del_gamma = del_q(length(k_tot));
    
    tol = 1e-8;
    tock = 0;
    error = abs(delta);
    
    while error > tol
        tock = tock + 1;
        u_bar = u_bar + delta_u;
        disp5 = [u_bar(1),u_bar(2),u_bar(3)];
        disp6 = [u_bar(4),u_bar(5),u_bar(6)];
        diff = disp6 - disp5;
        gamma = gamma + del_gamma;
        
        %%----Find Internal Force and stiffness at Updated Position-----%%
        F_int = [0,0,0,0,0,0];
        for m = 1:7 %m = element number
            if logic(m,3)==5
                [F{m},k{m}] = stiffness(EA, Length{m}, disp5, 3);
                F_int5 = F{m};
                F_int6 = [0,0,0];
            elseif logic(m,2)~=5
                [F{m},k{m}] = stiffness(EA, Length{m}, disp6, 3);
                F_int5 = [0,0,0];
                F_int6 = F{m};
            else
                [F{m},k{m}] = stiffness(EA, Length{m}, diff, 3);
                F_int5 = -F{m}; %points opposite of element n vector
                F_int6 = F{m};
            end
            F_int = F_int + [F_int5,F_int6];
            F36 = norm(F{4});
            F46 = norm(F{6});
        end
        kff1 = k{1} + k{2} + k{3} + k{7}; 
        kff2 = -k{7}; %Should use fixity values from layout instead of hardcoding
        kff3 = -k{7};
        kff4 = k{2} + k{4} + k{6} + k{7};
        k_tan = [kff1,kff2;kff3,kff4];
        k_tot = [k_tan,-Pref;-ek,0]; %Appended stiffness matrix
        
        %%-----Find Residuals-----%%
        R = -F_int' + Pref*gamma;
        g = ek*u_bar' - tick*delta;
        R_tilde = [R;g];
        error = norm(R_tilde);
        
        %%-----Calculate New Displacements-----%%
        del_q = k_tot\R_tilde;
        delta_u = [del_q(1),del_q(2),del_q(3),del_q(4),del_q(5),del_q(6)];
        del_gamma = del_q(length(k_tot));
        
        if tock > 1000
            break
        else
        end
    end
    u_disp5(tick) = u_bar(1);
    v_disp5(tick) = u_bar(2);
    u_disp6(tick) = u_bar(4);
    v_disp6(tick) = u_bar(5);
    w_disp6(tick) = u_bar(6);
    vert = vert + delta;
end
subplot(2,1,1)
grid on
hold on
plot(u_disp5,load_factor)
plot(v_disp5,load_factor)
plot(w_disp5,load_factor)
legend('U5_u','U5_v','U5_w')
xlabel('Displacement (m)')
ylabel('Load Factor, gamma')
title('Load-Displacement Curves for Free Nodes of Truss')
subplot(2,1,2)
grid on
hold on
plot(u_disp6,load_factor)
plot(v_disp6,load_factor)
plot(w_disp6,load_factor)
legend('U6_u','U6_v','U6_w')
title('Load-Displacement Curves for Free Nodes of Truss')
xlabel('Displacement (m)')
ylabel('Load Factor, gamma')

figure
plot(layout(5,2)+v_disp5,layout(5,3)+w_disp5,'o')
hold on
plot(layout(6,2)+v_disp6,layout(6,3)+w_disp6,'x')
legend('Node5','Node6')
title('y-z Planar Motion of Free Nodes')

figure
plot(layout(5,1)+u_disp5,layout(5,2)+v_disp5,'o')
hold on
plot(layout(6,1)+u_disp6,layout(6,2)+v_disp6,'x')
legend('Node5','Node6')
title('x-y Planar Motion of Free Nodes')


