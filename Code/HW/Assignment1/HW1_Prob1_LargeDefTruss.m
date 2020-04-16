%% CESG 506 HW1 - DIFFERENT FORMULATIONS FOR LARGE DEFORMATION PROBLEMS
clear;clc;
% PROBLEM SPECIFIC PARAMETERS
EA = 2100; %kN (Axial Stiffness)
W = 5.5; %m (Horizontal component of length, fixed)
H = 0.5; %m (Vertical component of length, varies)
L = sqrt(W^2+H^2); %m (Undeformed Length)
num_points = 1500;
u = linspace(0,10*H,num_points); %m (Change in vertical component)
l = sqrt(W^2+(H-u).^2); %m (Deformed Length)
N = [cos(atan(H/W));
     sin(atan(H/W))]; %Undeformed direction vector (Normalized)
n = [cos(atan((H-u)/W));
     sin(atan((H-u)/W))]; %Deformed direction vector (Normalized)

 
%%-----FOUR DIFFERENT STRAIN MEASURES-----%%
% ENGINEERING STRAIN
e = (l-L)/L;
% GREEN-LAGRANGE STRAIN
EE = 1/2*(l.^2-L^2)/L^2;
% HENKEY (NATURAL) STRAIN
eps = 1/2*log(l.^2/L^2);
% ALMANSI STRAIN
eps_a = 1/2*(l.^2-L^2)./l.^2;

%%-----INTERNAL FORCE MAGNITUDES-----%%
f_e = EA.*e;
f_EE = EA.*EE;
f_eps = EA.*eps;
f_eps_a = EA.*eps_a;

%%-----EQUILIBRIUM ON DEFORMED VS. UNDEFORMED CONFIGURATIONS-----%%
% UNDEFORMED CONFIGURATION
fu_e = f_e.*N;
fu_EE = f_EE.*N;
fu_eps = f_eps.*N;
fu_eps_a = f_eps_a.*N;

% DEFORMED CONFIGURATION
fd_e = f_e.*n;
fd_EE = f_EE.*n;
fd_eps = f_eps.*n;
fd_eps_a = f_eps_a.*n;

%%-----PLOTTING-----%%
% UNDEFORMED CONFIGURATION
Pu_e = zeros(1,num_points);
Pu_EE = zeros(1,num_points);
Pu_eps = zeros(1,num_points);
Pu_eps_a = zeros(1,num_points);

for i = 1:num_points 
    Pu_e(i) = fu_e(2,i);
    Pu_EE(i) = fu_EE(2,i);
    Pu_eps(i) = fu_eps(2,i);
    Pu_eps_a(i) = fu_eps_a(2,i);
end

hold on
plot(u,Pu_e,'--')
plot(u,Pu_EE,'--')
plot(u,Pu_eps,'--')
plot(u,Pu_eps_a,'--')

% DEFORMED CONFIGURATION
Pd_e = zeros(1,num_points);
Pd_EE = zeros(1,num_points);
Pd_eps = zeros(1,num_points);
Pd_eps_a = zeros(1,num_points);

for i = 1:num_points %Getting magnitudes of internal forces
    Pd_e(i) = fd_e(2,i);
    Pd_EE(i) = fd_EE(2,i);
    Pd_eps(i) = fd_eps(2,i);
    Pd_eps_a(i) = fd_eps_a(2,i);
end

plot(u,Pd_e)
plot(u,Pd_EE)
plot(u,Pd_eps)
plot(u,Pd_eps_a)

%LINEAR APPROXIMATION
k = EA/L*(N*N');
d{1} = [];
P = [];
for i = 1:num_points
    d{i} = [0 u(i)];
    P(i) = k(2,:)*d{i}';
end

plot(u,-P)
title('Exploring Nonlinearity')
xlabel('Displacement')
ylabel('Applied Force')
legend('','','','','deformed equilibrium1','deformed equilibrium2','deformed equilibrium3','deformed equilibrium4','Linear Approx')



