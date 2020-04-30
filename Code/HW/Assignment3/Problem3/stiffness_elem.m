function [F,k_e] = stiffness_elem(axial, Length, pos, dim)
%CALCULATES INDIVIDUAL ELEMENT STIFFNESS and FORCE VECTORS
%   Length is a (row) vector
%   axial = EA (axial stiffness for the element)
%   pos = displacement of free node from original position (row vector)
%   dim = dimension of element stiffness matrix (e.g. 3D)

len = Length + pos; %these are (row) vectors
l_mag = sqrt(dot(len,len)); %length magnitude
L_mag = sqrt(dot(Length,Length));
dir = len./l_mag; %unit direction (row) vector of deformed shape
stretch = l_mag/L_mag;
eps = log(stretch); %strain

f = axial*eps; %internal force magnitude
F_row = f.*dir;
F = F_row';
k_e = axial/l_mag.*(dir'*dir) + f/l_mag.*(eye(dim)-dir'*dir);



