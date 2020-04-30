function [k_append] = append(kff, Pref, u_bar, u_con, alpha, gamma, gamma_con)
%   APPENDS THE FREE-FREE SITFFNESS MATRIX WITH REFERENCE LOADS AND
%   ARC-LENGTH CONSTRAINT
%   

k_append = [kff, -Pref; 
        -2*(u_bar-u_con),-2*alpha*(gamma-gamma_con)*dot(Pref,Pref)];

end

