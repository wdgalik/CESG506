function u_full = rowfree2full(u_bar,idxF,idxS,numDOF,numNodes)

%   INPUTS
%   u_bar: column vector of free displacements
%   idxF:
%   idxS:

%   OUTPUTS
%   u_full: displacement array (numDOF x numNodes)


u_col = zeros(numDOF*numNodes,1);
for i = 1:length(idxS)
    u_col(idxS(i)) = 0;
end

for j = 1:length(idxF)
    u_col(idxF(j)) = u_bar(j);
end

u_full = reshape(u_col,numDOF,numNodes);
end

