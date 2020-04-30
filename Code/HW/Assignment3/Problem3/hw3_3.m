%% CESG 506 HW2 - ARC LENGTH METHOD - TRUSS COLUMN
clear; clc;
%%-----PROBLEM SPECIFIC PARAMETERS-----%%
EA_vert = 2000; %kN
EA_brace = 5000;
DOF = 2; %# number nodal dof

nodal_loads = [23, -3, 2;
               24, -3, 2]; %(node#, value, dof)

%%-----DEFINE SYSTEM LAYOUT-----%%
num_bays = 11;
num_nodes = (num_bays+1)*2;
num_per_bay = 4;
num_elem = num_bays*num_per_bay;
H = 5;
H_bay = H/num_bays;
W = H/20;

%%-----Define Nodal Coordinates-----%%
layout = cell(num_nodes,1);
minus = 1;
sub = 2;
for mm = 1:num_nodes
    if mod(mm, DOF) == 1
        layout{mm} = [0, (mm-minus)*H_bay,1]; %x,y,fixity (0=fixed, 1=free)
        minus = minus + 1;
    elseif mod(mm, DOF) == 0
        layout{mm} = [W, (mm-minus)*H_bay,1];
        sub = sub + 1;
    end
end
layout{1}(DOF+1) = 0;
layout{2}(DOF+1) = 0;
coordinates = cell2mat(layout); %nodal layout


%%-----Define Element Connectivity-----%%
mesh = cell(num_elem,1); 
subi = [1, 2, 1, 3];
subj = [3, 4, 4, 4];
for ii = 1:num_bays
    less = 3;
    for jj = 1:num_per_bay
        if jj == 1 || jj == 2
            index = ii*num_per_bay-less;
            mesh{index} = [index, subi(jj), subj(jj), EA_vert];
            less = less - 1;
        else
            index = ii*num_per_bay-less;
            mesh{index} = [index, subi(jj), subj(jj), EA_brace];
            less = less - 1;
        end
    end
    subi = subi + 2;
    subj = subj + 2;
end
logic = cell2mat(mesh); %Element connectivity (element#,nodei, nodej, EA)

     
for i = 1:length(logic)
    Length{i} = coordinates(logic(i,3),(1:DOF))-coordinates(logic(i,2),(1:DOF));
end

%%-----INITIALIZATION-----%%
[disp, u_bar, num_free, node_free, R, F_int] = init(coordinates, DOF);

Pcr = pi^2*EA_vert*W^2/8/H^2;
Pref = zeros(num_free*DOF, 1);
Pref(num_free*DOF) = -Pcr/2;
Pref(num_free*DOF-2) = -Pcr/2;

u_con = u_bar; %initial converged displacement
u_storage = {};
u_storage = [u_storage, u_con];
gamma_con = 0;
gamma = gamma_con; %load factor
gamma_storage = [];
gamma_storage = [gamma_storage, gamma_con];


%%-----Get Initial Element Stiffnesses-----%%
ss = size(logic);
for m = 1:ss(1) %m = element number
    diff = cell2mat(disp(logic(m,3)))-cell2mat(disp(logic(m,2)));
    EA = logic(m,DOF+2);
    [F{m},k{m}] = stiffness_elem(EA, Length{m}, diff, DOF);
end

%%-----Assembling/Partitioning Full Stiffness Matrix
[kff] = assembly(k, coordinates, logic, num_free, node_free, DOF);


%%----ITERATIVE NEWTON-RAPHSON ALGORITHM (ARC LENGTH METHOD)-----%%

%%-----Determine Step size, delta_s-----%%
alpha = 0.05;
dgamma = 0.2;
du = kff\(Pref*dgamma);
delta_s = sqrt(dot(du,du)+alpha*dgamma^2*dot(Pref,Pref));

tick = 0;
while gamma < 2
    tick = tick + 1;
    %%-----Converged Conditions-----%%
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
    
    tol = 0.00001;
    tock = 0;
    error = delta_s^2; 
    while error > tol
        tock = tock + 1;
        
        for p = 1:num_free
           st = DOF*p-(DOF-1);
           disp{node_free(p)} = u_bar(st:DOF*p);
        end
        
        %%----Find Internal Force and stiffness at Updated Position-----%%
        for ii = 1:num_free
            F_int{ii} = zeros(DOF,1);
        end
        for n = 1:length(logic(:,1)) %n = element number
            diff = cell2mat(disp(logic(n,3)))-cell2mat(disp(logic(n,2)));
            EA = logic(n,DOF+2);
            [F{n},k{n}] = stiffness_elem(EA, Length{n}, diff, DOF); 
            for b = 1:num_free
               if logic(n,3) == node_free(b)
                   F_int{b} = F_int{b} + F{n};
               elseif logic(n,2) == node_free(b)
                   F_int{b} = F_int{b} - F{n};
               else
               end
            end
        end
        
        %%-----Assembling/Partitioning Updated Stiffness Matrix-----%%
        [kff] = assembly(k, coordinates, logic, num_free, node_free, DOF);
        
        %%-----Appending Stiffness Matrix-----%%
        [k_tot] = append(kff, Pref, u_bar, u_con, alpha, gamma, gamma_con);
        
        %%-----Find Residuals-----%%
        R = -cell2mat(F_int') + Pref*gamma;
        g = 0;
        R_tilde = [R;g];
        error = norm(R_tilde);
        
        %%-----Calculate New Displacements-----%%
        del_q = k_tot\R_tilde;
        u_bar = u_bar + del_q(1:num_free*DOF)';
        gamma = gamma + del_q(length(k_tot));
        if tock > 100
            break
        else
        end
    end
    top2y_disp(tick) = u_storage{tick}(num_free*DOF);
    top2x_disp(tick) = u_storage{tick}(num_free*DOF-1);
    top1y_disp(tick) = u_storage{tick}((num_free-1)*DOF);
    top1x_disp(tick) = u_storage{tick}((num_free-1)*DOF-1);
end

%%-----Plotting Load-Deformation Curves for top two nodes-----%%
plot(top2x_disp, gamma_storage,'--')
hold on
plot(top1x_disp, gamma_storage,'--')
plot(top2y_disp, gamma_storage)
plot(top1y_disp, gamma_storage)
xlabel('Displacement, m')
ylabel('Applied Force, kN')
title('Load-Displacement Curves')
legend('x-disp: No.24','x-disp: No.23','y-disp: No.24','y-disp: No.23')


%%-----Path Tracing for Top Two Nodes-----%%
figure
plot(top1x_disp,top1y_disp + coordinates(24,2),'--')
hold on
plot(top2x_disp + coordinates(24,1),top2y_disp + coordinates(24,2),'--')
xlabel('x-displacement')
ylabel('y-displacement')
title('Deformed System')

for nn = 1:num_nodes
    x_rest(nn) = disp{nn}(1) + coordinates(nn,1);
    y_rest(nn) = disp{nn}(2) + coordinates(nn,2);
end

for qq = 1:num_elem    
    plot([x_rest(logic(qq,2)) x_rest(logic(qq,3))],[y_rest(logic(qq,2)) y_rest(logic(qq,3))])
    hold on
end
   

    
