function [Pf, Ff, Kff] = beam_assembly(nodes, mesh, fixity, P, U)

%   INPUTS
%   nodes: matrix that describes nodal coords (node ID, x, y)
%   mesh: element connectivity (element#, i, j, EA, EI)
%   fixity: describes support state (node#, x, y, theta) 1=fixed, 0=free
%   P: Applied force (needs to be full size...not reduced to free dofs)
%   U: nodal displacements [u1 u2...un; v1 v2...vn; theta1,theta2...thetan]
        ....needs to be full size (numDOF x numNodes)

%   OUTPUTS
%   Psys: Applied force vector
%   Fsys: Internal force vector
%   Ksys: System stiffness matrix

numDOF = 3; %(u,v,theta)
numNodes = length(nodes(:,1));
fullsize = numDOF*numNodes;

% Initialize internal force vector and global stiffness matrix
%   Define cells at each node to help with assembly
F_rowsplit = numDOF*ones(1,numNodes);
FsysC = mat2cell(zeros(fullsize,1),[F_rowsplit],1); 

K_rowsplit = numDOF*ones(1,numNodes);
K_colsplit = numDOF*ones(1,numNodes);
KsysC = mat2cell(zeros(numDOF*numNodes),K_rowsplit,K_colsplit);

% get initial nodal coordinates
ncoords = nodes(:,2:end)'; %[x1 x2 x3...xn; y1 y2 y3...yn];

% loop through all elements and 'grab' properties
for elem = 1:length(mesh(:,1))
	i  = mesh(elem,2);
	j  = mesh(elem,3);
	EA = mesh(elem,4);
	EI = mesh(elem,5);

    % get fe, ke
	[fe,ke] = curved_beam_elem(ncoords(:,i), ncoords(:,j), U(:,i), U(:,j), EA, EI);

    % assemble element e to system (C for cell)
	FsysC{i} = FsysC{i} + fe(1:3);
	FsysC{j} = FsysC{j} + fe(4:6);

	KsysC{i,i} = KsysC{i,i} + ke(1:3, 1:3);
	KsysC{i,j} = KsysC{i,j} + ke(1:3, 4:6);
	KsysC{j,i} = KsysC{j,i} + ke(4:6, 1:3);
	KsysC{j,j} = KsysC{j,j} + ke(4:6, 4:6);
end % end element loop

% cell to matrix (for practicality)
Psys = reshape(P',fullsize,1);
Fsys=cell2mat(FsysC);
Ksys=cell2mat(KsysC);

% reduce problem to free dofs
[Pf, Ff, Kff] = reduction(numDOF, numNodes, fixity, Psys, Fsys, Ksys);


end % end function definition

