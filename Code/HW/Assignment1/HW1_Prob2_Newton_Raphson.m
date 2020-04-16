%% CESG 506 HW1 - NEWTON RAPHSON ITERATION METHOD
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
Pcr = [0; 0.9817]; %kN 
gamma = [0,0.25,0.5,0.75,0.99,0.999]; %load factors for iterative NR method

%%-----ITERATIVE NEWTON RAPHSON ALGORITHM-----%%
Residual = [];
for i = 1:(length(gamma)-1)
    if i == 1 %starting iteration at undeformed configuration
        u_bar = [0,0];
        ticker = 0;        
        R = gamma(2)*Pcr; %for zero displacement, internal force = 0
        Residual(ticker+1,i) = norm(R);
        k_init = EA/L1*(L1_vec'*L1_vec) + EA/L2*(L2_vec'*L2_vec); %note that 2nd term absent b/c int force = 0 for no init disp
        h = -inv(k_init)*R;
    else
        ticker = 0;
        k_init = k_tan; %k_tan calculated in first while loop
        R = gamma(i+1)*Pcr + F_int';
        Residual(ticker+1,i) = norm(R);
        h = -inv(k_init)*R;
    end
tol = 1e-15;
while norm(R) > tol
    ticker = ticker + 1;
    %%-----UPDATE POSTION VECTOR-----%%
    u_bar = u_bar + h'; 
    
    %%----FIND INTERNAL FORCE VECTOR AT UPDATED POSITION-----%%
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
    
    %%-----CALCULATE RESIDUAL FORCE FOR UPDATED POSITION-----%%
    R = gamma(i+1)*Pcr + F_int';
    Residual(ticker+1,i) = norm(R);
    error = norm(R);
    
    %%-----CALCULATE UPDATED TANGENT STIFFNESS-----%%
    k1 = (EA/l1_mag).*(n_1'*n_1) + (f_1/l1_mag).*(eye(2)-n_1'*n_1);
    k2 = EA/l2_mag.*(n_2'*n_2) + f_2/l2_mag.*(eye(2)-n_2'*n_2); 
    k_tan = k1 + k2;
    
    %%-----CALCULATE NEW DISPLACEMENT DELTA (h)-----%%
    h = -inv(k_tan)*R;
    
    if ticker > 25
        break
    else
    end
end
% formatSpec = 'Gamma is %3.2f \n';
% fprintf(formatSpec,gamma)
line = '\n Displacement ';
number = 'at load level %1.0f ';
is = 'is: \n';
print = '%6.5f \n';
fprintf(line)
fprintf(number,i+1)
fprintf(is)
fprintf(print,u_bar)
end

fprintf('\n TABLE OF RESIDUALS: \n')
T = array2table(Residual);
disp(T)




    
    
    
        
        