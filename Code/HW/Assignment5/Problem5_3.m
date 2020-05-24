%% Problem 5-3 SHALLOW ARCH SNAP-THROUGH
clear;clc;

%%----------------------%%
%%---PROBLEM SPECIFIC---%%
%%----------------------%%
W = 5; %total width (m)
H = W/10; %midspan height (m)

w = 1; %distributed load (kN/m)

EA  = 10000; %kN
EI  = 10; %kN-m^2


%%-----------------------%%
%%----DEFINE MESHING-----%%
%%-----------------------%%
numElems = 64;
nNds = numElems + 1;

nDOF = 3; %u,v,theta

nodeID = linspace(1, nNds, nNds);
elemID = linspace(1, numElems, numElems);
X  = linspace(0, W, nNds); % define x position of each node
x_norm = X/W; % normalized x-position
h  = 4*H.*x_norm.*(1-x_norm); % y-position of each node

nodes = [nodeID', X', h' ]; %nodeID, x, y

mesh = [elemID', elemID', elemID'+1, EA*ones(numElems,1), EI*ones(numElems,1) ];
        
fixity = [1    1    1   0; % node, x, y, theta
         nNds  1    1   0];

%%---------------------------%%
%%-----APPLY NODAL LOADS-----%%
%%---------------------------%%
Le = W/numElems;
Pe = w*Le/2; %force at each node of element due to distr. load
Me = w*Le*Le/12; 

loads = [nodeID', zeros(nNds, nDOF)]; % node, Px, Py, M
loads(:,3)    = -2*Pe;  % interior nodes see forces from 2 elements
loads(1,3)    =   -Pe;  % first node sees force from first element only
loads(nNds,3) =   -Pe;  % last node sees force from last element only
loads(1,4)    =   -Me;  % moments cancel on interior nodes but not at either end
loads(nNds,4) =    Me;  % moments cancel on interior nodes but not at either end

% create the imperfection
lf = 0.99; %relative magnitude of right side dist. load
rf = 1.01; %relative magnitude of left side dist. load
factors = 1.0 + (rf-lf)*((nodeID - 1)/(nNds - 1) - 1/2 );
loads(:,3) = loads(:,3) .* factors';

Pfullsize = loads(:,2:end);

%%------------------------%%
%%-----INITIALIZATION-----%%
%%------------------------%%
U0 = zeros(nDOF,nNds); %initial (full) displacement matrix of zeros
[Pref, F_init, K_init] = beam_assembly(nodes, mesh, fixity, loads(:,2:end), U0);

% create storage cells
full_disp = reshape(U0,nDOF*nNds,1); %initial full displ vector
u_bar = full_disp(1:length(Pref),:); %initial free displ vector
u_con = u_bar; %initial converged displacement
u_storage = {};
u_storage = [u_storage, u_con];
gamma_con = 0;
gamma = gamma_con; %load factor
gamma_storage = [];
gamma_storage = [gamma_storage, gamma_con];

%%-----------------------------------------------------%%
%%-----NEWTON RAPHSON ALGORITHM - ARCLENGTH METHOD-----%%
%%-----------------------------------------------------%%
% Determine Step size, delta_s
alpha = 0.0;
dgamma = 0.2;
du = K_init\(Pref*dgamma);
delta_s = sqrt(dot(du,du)+alpha*dgamma^2*dot(Pref,Pref));

tick = 0;
while gamma >= 0
    tick = tick + 1;
    %%-----Converged Conditions-----%%
    if tick == 1
        u_bar = u_bar + du;
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
        error = delta_s; 
        
        while error > tol
            tock = tock + 1;
            
        %%----Find Internal Force and stiffness at Updated Position-----%%
        [idxF, idxS] = get_indices(nDOF, nNds, fixity);
        u_full = rowfree2full(u_bar,idxF,idxS,nDOF,nNds);
        [P, F_int, kff] = beam_assembly(nodes, mesh, fixity, Pfullsize*gamma, u_full);  
        
        %%-----Appending Stiffness Matrix-----%
        [k_tot] = append(kff, Pref, u_bar', u_con', alpha, gamma, gamma_con);
        
        %%-----Find Residuals-----%%
        R = -F_int + Pref*gamma;
        g = 0;
        R_tilde = [R;g];
        error = norm(R_tilde);
        
        %%-----Calculate New Displacements-----%%
        del_q = k_tot\R_tilde;
        u_bar = u_bar + del_q(1:length(kff));
        gamma = gamma + del_q(length(k_tot));
        
        if tock > 25
            break
        else
        end
        
        end %end inner loop of iterations
        
end %end outer loop of steps

%%------------------------------------------%%
%%-----PLOTTING DISPLACED SHAPE RESULTS-----%%
%%------------------------------------------%%
% scale free displacements back up to full size
u_final = zeros(nDOF*nNds,1);

% assign fixed dofs to have zero displacement
for s = 1:length(idxS)
    u_final(idxS(s)) = 0;
end

% parse in the free displacements into full vector
for f = 1:length(idxF)
    u_final(idxF(f)) = u_con(f);
end

% extract vectors for u and v displacements
del_x_nodes = u_final(1:3:end);
del_y_nodes = u_final(2:3:end);

% final position of the nodes
x_end = X + del_x_nodes';
y_end = h + del_y_nodes';

% overlay deformed and undeformed configurations
plot(X,h)
hold on
plot(x_end,y_end)
title('Beam Plot')
legend('Undeformed','Deformed')
xlabel('x (m)')
ylabel('h (m)')

%%------------------------------------------------------%%
%%-----PLOTTING LOAD FACTORS VS NODAL DISPLACEMENTS-----%%
%%------------------------------------------------------%%
% get full displacement vector for each converged step
u_store_full = cell(1,tick);
for it = 1:tick
    store_full = zeros(nDOF*nNds,1);
    for sup = 1:length(idxS)
        store_full(idxS(sup)) = 0;
    end  
    for fr = 1:length(idxF)
        store_full(idxF(fr)) = u_storage{it}(fr);
    end
    u_store_full{it} = store_full;
end

% extract nodal displacements from u_storage cells
mid_node = (nNds+1)/2;
quarter_node = ((nNds+1)/2+1)/2;
q3_node = mid_node + quarter_node;
v_mid = zeros(1,tick);
v_quarter = zeros(1,tick);
v_q3 = zeros(1,tick);
for t = 1:tick
    v_mid(t) = u_store_full{t}(nDOF*(mid_node-1)+2);
    v_quarter(t) = u_store_full{t}(nDOF*(quarter_node-1)+2);
    v_q3(t) = u_store_full{t}(nDOF*(q3_node-1)+2);
end

% compare y-displacements of middle and quarter points
figure
plot(v_mid,gamma_storage)
hold on
plot(v_quarter,gamma_storage)
plot(v_q3,gamma_storage)
grid on
legend('U_v 1/4 Node','U_v 1/2 Node','U_v 3/4 Node')
title('Load Factor vs. Nodal y-displacements')
xlabel('Displacement (m)')
ylabel('Load Factor')




