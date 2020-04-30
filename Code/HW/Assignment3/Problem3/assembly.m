function [k_ff] = assembly(k_e, layout, connection, num_free, node_free, dof)
%BUILDS FULL STIFFNESS AND PARTITIONS
%   Input each element stiffness matrix, the nodal layout and connection
%   scheme
%   Output kff and appended stiffness matrix

%%-----Assembling Full Stiffness Matrix-----%%
s = size(layout);
k_full = cell(s(1),s(1));
for w = 1:s(1)
    for y = 1:s(1)
        k_full{w,y} = zeros(dof);
    end
end

for r = 1:length(connection)
    k_full{connection(r,2),connection(r,2)} = k_full{connection(r,2),connection(r,2)} + k_e{r};
    k_full{connection(r,2),connection(r,3)} = k_full{connection(r,2),connection(r,3)} - k_e{r};
    k_full{connection(r,3),connection(r,2)} = k_full{connection(r,3),connection(r,2)} - k_e{r};
    k_full{connection(r,3),connection(r,3)} = k_full{connection(r,3),connection(r,3)} + k_e{r};
end

%%-----Partitioning to get Kff-----%%
kff = [];
for q = 1:num_free
    for p = 1:num_free
        kff{q,p} = k_full{node_free(q),node_free(p)};
    end
end
k_ff = cell2mat(kff);

