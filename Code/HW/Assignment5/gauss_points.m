function [eval_pos,eval_weight] = gauss_points(ngp)

%   INPUTS
%   ngp: number of gauss points (evaluate polynomial of O(2n-1)
%
%   OUTPUTS
%   eval_pos: gives array containing gauss integration evaluation locations
%   eval_weight: corresponding weight factors for gauss integration

%   ***NOTE**** typical gauss integration tables are made for the interval
%   from [-1,1]. The beam element has length normalized from [0,1]. This is
%   reflected in the following evaluation locations. 

%   ***NOTE*** this function is limited to 5 gauss point integration
%   accuracy

%converting weights from normal gaussian quadrature [-1,1] to altered range
%[0,1]
Le = 1;
typ_range = 2; 
wf = Le/typ_range;

if ngp == 1 
    eval_pos = 0.5;
    eval_weight = 2/wf;
    
elseif ngp == 2
    typ_gp1 = -1/sqrt(3);
    typ_gp2 = 1/sqrt(3);
    
    eval_pos = [Le/2*(1 + typ_gp1), Le/2*(1 + typ_gp2)];
    eval_weight = [1,1]./wf;
    
elseif ngp == 3
    typ_gp1 = -3/sqrt(5);
    typ_gp2 = 0;
    typ_gp3 = 3/sqrt(5);
    
    eval_pos = [Le/2*(1 + typ_gp1),Le/2*(1 + typ_gp2),Le/2*(1 + typ_gp3)];
    eval_weight = [5/9 8/9 5/9]./wf;
    
elseif ngp == 4
    typ_gp1 = -sqrt((3+2*sqrt(6/5))/7);
    typ_gp2 = -sqrt((3-2*sqrt(6/5))/7);
    typ_gp3 = sqrt((3-2*sqrt(6/5))/7);
    typ_gp4 = sqrt((3+2*sqrt(6/5))/7);
    
    eval_pos = [Le/2*(1 + typ_gp1),Le/2*(1 + typ_gp2),Le/2*(1 + typ_gp3),Le/2*(1 + typ_gp4)];
    eval_weight = [(18-sqrt(30))/36,(18+sqrt(30))/36,(18+sqrt(30))/36,(18-sqrt(30))/36]./wf;
    
else 
    typ_gp1 = -sqrt(5+2*sqrt(10/7))/3;
    typ_gp2 = -sqrt(5-2*sqrt(10/7))/3;
    typ_gp3 = 0;
    typ_gp4 = sqrt(5-2*sqrt(10/7))/3;
    typ_gp5 = sqrt(5+2*sqrt(10/7))/3;
    
    eval_pos = [Le/2*(1 + typ_gp1),Le/2*(1 + typ_gp2),Le/2*(1 + typ_gp3),Le/2*(1 + typ_gp4),Le/2*(1 + typ_gp5)];
    eval_weight = [(322-13*sqrt(70))/900,(322+13*sqrt(70))/900,128/225,(322+13*sqrt(70))/900,(322-13*sqrt(70))/900./wf
    
end

