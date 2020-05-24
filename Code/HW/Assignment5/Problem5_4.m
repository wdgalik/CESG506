%% PROBLEM 5-4 BEAM BUCKLING
clear;clc;

%%--------------------------%%
%%-----PROBLEM SPECIFIC-----%%
%%--------------------------%%
% geometry
L = 5.00;
% H = L/100; % imperfect beam
H = 0; % straight beam

% stiffness properties
EA  = 10000; 
EI  = 10;

%%-----------------------%%
%%----DEFINE MESHING-----%%
%%-----------------------%%
numElems = 4;

nDOF = 3;
nNds = numElems + 1;

nodeID = linspace(1, nNds, nNds);
elemID = linspace(1, numElems, numElems);
X  = linspace(0, L, nNds);
x_norm = X/L;
h  = 4*H.*x_norm.*(1-x_norm);

% define initial nodal layout (nodeID, x, y)
nodes = [ nodeID', X', h' ];

% define element connectivity (elemID, i, j, EA, EI)
mesh = [  elemID', elemID', elemID'+1, EA*ones(numElems,1), EI*ones(numElems,1)];
        
% define support condition (node, x, y, theta) 1=fix, 0=free
fixity = [  1       1       1      0; 
           nNds     0       1      0]; 

%%---------------------------%%
%%-----APPLY NODAL LOADS-----%%
%%---------------------------%%
P = pi^2*EI/(L^2); %euler buckling load

loads=[nodeID', zeros(nNds, nDOF)]; %node, Px, Py, M
loads(nNds,2) =  -P ;  % apply axial load at the end

Pfullsize = loads(:,2:end);

%%------------------------%%
%%-----INITIALIZATION-----%%
%%------------------------%%
U0 = zeros(nDOF,nNds); %initial (full) displacement matrix of zeros
[Pref, F_init, K_init] = beam_assembly(nodes, mesh, fixity, Pfullsize, U0);

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
% Determine Step size, delta_s - constant throughout algorithm
alpha = 0.0;
dgamma = 0.01;
du = K_init\(Pref*dgamma);
delta_s = sqrt(dot(du,du)+alpha*dgamma^2*dot(Pref,Pref));

% begin outer loop - converged steps
tick = 0;
while gamma < 1.5
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

    % begin inner loop - iterate for convergence
        tol = 0.00001;
        tock = 0;
        error = delta_s; 
        while error > tol
            tock = tock + 1;
            
        %%----Find Internal Force and stiffness at Updated Position-----%%
        [idxF, idxS] = get_indices(nDOF, nNds, fixity);
        u_full = rowfree2full(u_bar,idxF,idxS,nDOF,nNds);
        [P, F_int, kff] = beam_assembly(nodes, mesh, fixity, Pfullsize*gamma, u_full);
        
        %%-----Get det of tangent stiffness (check buckled)-----%%
        detk(tick) = det(kff);
        
        %%-----Appending Stiffness Matrix-----%
        k_tot = append(kff, Pref, u_bar', u_con', alpha, gamma, gamma_con);
        
        %%-----Find Residuals-----%%
        R = -F_int + Pref*gamma;
        g = 0; %arc length constraint is satisfied for any displ guess
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

%%----------------------%%
%%-----PLOT RESULTS-----%%
%%----------------------%%
horiz = zeros(tick,1);
plot(gamma_storage, detk)
hold on
plot(gamma_storage, horiz)
