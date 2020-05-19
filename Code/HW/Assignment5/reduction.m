function [Pf, Ff, Kff] = reduction(nDOF, nds, fixity, Psys, Fsys, Ksys)

%   INPUTS
%   nDOF: number of dofs per node
%   nds: number of nodes
%   fixity: matrix that describes support conditions at nodes (free/fixed)
%   Psys: full external force array
%   Fsys: full internal force array
%   Ksys: full tangent stiffness matrix

%   OUTPUTS
%   Pff: Applied forces at free nodes
%   Fff: Intenal forces at free nodes
%   Kff: tangent stiffness matrix partitioned at free dofs

% intially, set all nodes as free
idxF = linspace(1,nDOF*nds,nDOF*nds);
idxS = [];

% find supported nodes
for nsn = 1:length(fixity(:, 1))
	nodeID = support(nsn, 1);
    for dof = 1:nDOF
	    if support(nsn, dof+1) ~= 0
            idxS(end+1) = (nodeID-1)*nDOF + dof;
	    end
    end
end

% remove supported nodes from free list
idxF(idxS) = [];

% reduce the stiffness matrix and the residual to free DOFs only
Pf = Psys;
Ff = Fsys;
Kff = Ksys;

Pf(idxS)  =[];
Ff(idxS)  =[];
Kff(:,idxS)=[];
Kff(idxS,:)=[];

end

