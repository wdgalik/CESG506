function shapes = shape_function(order,d,x_norm,L_elem)

%   INPUTS
%   order: tells order of shape function
%   d: how many derivatives to take of shape function
%   x_norm: normalized position along length of element
%   L_elem: Length of element
%   
%   OUTPUTS
%   shapes: vector of shape functions

if order == 1 %linear displacement (u) shape functions
    if d == 0
        shapes = [1-x_norm, x_norm];
    elseif d ==1
        shapes = [-1/L_elem, 1/L_elem];
    else
        shapes = [0, 0];
    end
elseif order == 3 %cubic displacement (v) shape functions
    if d == 0
        shapes = [1-3*x_norm^2+2*x_norm^3, L_elem*(x_norm-2*x_norm^2+x_norm^3),...
                  3*x_norm^2-2*x_norm^3, L_elem*(x_norm^3-x_norm^2)];
    elseif d ==1
        shapes = [-6*x_norm+6*x_norm^2, L_elem*(1-4*x_norm+3*x_norm^2),...
                  6*x_norm-6*x_norm^2, L_elem*(-2*x_norm+3*x_norm^2)]/L_elem;
    elseif d ==2
        shapes = [-6+12*x_norm, L_elem*(-4+6*x_norm), 6-12*x_norm,...
                  L_elem*(-2+6*x_norm)]/(L_elem^2);
    elseif d == 3
        shapes = [12, L_elem*6, -12, L_elem*6]/(L_elem^3);
    else
        shapes = [0,0,0,0];
    end
end

