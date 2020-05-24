function [idxF, idxS] = get_indices(nDOF, nNds, fixity)

%   INPUTS
%   nDOF: number of dofs per node
%   nds: number of nodes
%   fixity: matrix that describes support conditions at nodes (free/fixed)

%   OUTPUTS
%   idxF: indices of free nodes
%   idxS: indices of supported nodes

    % initialize all as free
    idxF = linspace(1,nDOF*nNds,nDOF*nNds);
    idxS = [];

    % find supported node indices (nsn = num supported nodes)
    for nsn=1:length(fixity(:,1))
	nodeID = fixity(nsn,1);
	for i=1:nDOF
	    if fixity(nsn,i+1)>0
		idxS(end+1) = nodeID*nDOF-nDOF+i;
	    end
	end
    end

    % remove supported nodes from free list
    idxF(idxS) = [];

end