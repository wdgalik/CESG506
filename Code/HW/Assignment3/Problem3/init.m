function [disp_c, u_free, num_free, node_free, R, F_int] = init(layout, dof)
%TAKES NODAL LAYOUT AND INITIALIZES DISPLACEMENTS, RETURNS NUMBER OF FREE
%NODES AND RETURNS FREE NODAL NUMBERS
%   note that this function assumes no support settlements
%   all displacements initialized to zero

s = size(layout);
num_free = 0;
node_free = [];
u_free = [];

for i = 1:s(1)
    disp_c{i} = zeros(1,s(2)-1); %initializing nodal displacements
    if layout(i,s(2)) == 1 %checking fixity
        num_free = num_free + 1;
        node_free = [node_free, i];
        u_free = [u_free, disp_c{i}];
    else 
    end
end

free_dof = cell(num_free*dof,1);
R = zeros(num_free*dof,1);
F_int = cell(1,num_free);
for ii = 1:num_free
    F_int{ii} = zeros(1,dof);
end


