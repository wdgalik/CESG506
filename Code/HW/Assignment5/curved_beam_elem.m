function [Fe, Ke, dNv] = curved_beam_elem(Pos_i, Pos_j, dof_i, dof_j, EA, EI)

%   INPUTS
%   Pos_i: position vector at node i (xi, hi)
%   Pos_j: position vector at node j (xj, hj)
%   dof_i: dofs (displacement values) at node i (ui, vi, thetai)
%   dof_j: dofs (displacement values) at node j (uj, vj, thetaj
%   EA: Axial stiffness
%   EI: Bending Stiffness

%   OUTPUTS
%   Fe: Force vector (Internal)
%   Ke: Tangent stiffness matrix

Le = Pos_j(1) - Pos_i(1);

% Partition nodal dofs
qu = [dof_i(1); dof_j(1)];
qv = [dof_i(2); dof_i(3); dof_j(2); dof_j(3)];

% Initialize force and stiffness
Fe = zeros(6,1);
Ke = zeros(6,6);
% P = zeros(6,1);

% Gaussian-integration parameters
ngp = 2;
[eval_pos,eval_weight] = gauss_points(ngp);

% Start numerical integration loop
for gp = 1:ngp
    % get function evaluation location and weighting value
    posit = eval_pos(gp);
    weight = eval_weight(gp)*Le;
    
    % obtain shape functions
    Nu = shape_function(1,0,posit,Le);
    dNu = shape_function(1,1,posit,Le);
    Nv = shape_function(3,0,posit,Le);
    dNv = shape_function(3,1,posit,Le);
    ddNv = shape_function(3,2,posit,Le);
    
    % kinematics (midline approximations)
    dh = (Pos_j(2) - Pos_i(2))/Le;
    du = dNu*qu;
    dv = dNv*qv;
    phi = ddNv*qv;
    eps = du + dh*dv + 0.5*dv*dv;
    
    % internal resultants
    f = EA*eps;
    M = EI*phi;
    
    % Strain-Displacement, B, Matrices (provide a linear map between..
    %    ..variations of nodal dofs and variation of strains
    B_eps = [dNu(1), (dh+dv)*dNv(1), (dh+dv)*dNv(2), dNu(2),...
            (dh+dv)*dNv(3), (dh+dv)*dNv(4)];
    B_phi = [0, ddNv(1), ddNv(2), 0, ddNv(3), ddNv(4)];
    
    % Define Nv_hat, which is like Nv, but is in same dimension as dofs
    Nv_hat = [0, Nv(1), Nv(2), 0, Nv(3), Nv(4)];
    dNv_hat = [0, dNv(1), dNv(2), 0, dNv(3), dNv(4)];
    
    % Get internal force vector
    Fe = Fe + f*B_eps'*weight + M*B_phi'*weight;
    
    % Get applied force vector
%     P = P - w.*Nv_hat'*weight;

    % Get element stiffness matrix
    Ke = Ke + EA*(B_eps'*B_eps)*weight + EI*(B_phi'*B_phi)*weight + ...
         f*(dNv_hat'*dNv_hat)*weight; 
    
end %end gauss integration loop

end

