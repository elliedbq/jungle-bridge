% day 21

%computes the distance constraint error for a pair of vertices
    %INPUTS:
        %(xA, yA): coordinates of first vertex
        %(xB, yB): coordinates of second vertex
        %l_max: maximum allowable distance between two vertices
    %OUTPUTS:
        %e_len: constraint error.
        % e_len<=0 when constraint is satisfied
        % e_len>0 when constraint is violated
function e_len = single_string_error_func(xA,yA,xB,yB,l_max)
    e_len = sqrt((xB-xA)^2 + (yB-yA)^2) - l_max;
end


%evaluates the distance constraint error across all links
%INPUTS:
%coords: vector of vertex positions from i=1 to i=(n-1)
% [x_1;y_1;...;x_(n-1),y_(n-1)]
%param_struct: struct containing parameters of the bridge
% param_struct.r0 = [x_0;y_0]: coordinates of leftmost vertex
% param_struct.rn = [x_n;y_n]: coordinates of rightmost vertex
% param_struct.num_links: number of links in bridge
% param_struct.l0_list = [l0_1;...;l0_n]: list of link lengths
% param_struct.m_list = [m_1;...;m_(n-1)]: list of weight masses
% param_struct.g = 9.8 m/sec^2: gravitational acceleration
%OUTPUTS:
%e_val = [e_len1; ... ; e_len_n]: the vector of distance constraint errors
%dummy = []: empty vector used to satisfy fmincon syntax
%evaluates the distance constraint error across all links
%INPUTS:
%coords: vector of vertex positions from i=1 to i=(n-1)
% [x_1;y_1;...;x_(n-1),y_(n-1)]
%param_struct: struct containing parameters of the bridge
% param_struct.r0 = [x_0;y_0]: coordinates of leftmost vertex
% param_struct.rn = [x_n;y_n]: coordinates of rightmost vertex
% param_struct.num_links: number of links in bridge
% param_struct.l0_list = [l_1;...;l_n]: list of link lengths
% param_struct.m_list = [m_1;...;m_(n-1)]: list of weight masses
% param_struct.g = 9.8 m/sec^2: gravitational acceleration
%OUTPUTS:
%e_val = [e_len1; ... ; e_len_n]: the vector of distance constraint errors
%dummy = []: empty vector used to satisfy fmincon syntax
function [e_vec,dummy] = bridge_error_func(coords,param_struct)
    %initialize error vector
    e_vec = zeros(param_struct.num_links,1);
    %initialize dummy output for fmincon equality constraints
    dummy = [];
    %add the first and last vertex positions to the coordinate list
    coords = [param_struct.r0;coords;param_struct.rn];
    %iterate through each rubber band link
    for i = 1:param_struct.num_links
        %extract the ith segment length
        l_max = param_struct.l0_list(i);
        %extract the coordinates of the string ends
        xA = coords(2*i-1);
        yA = coords(2*i);
        xB = coords(2*i+1);
        yB = coords(2*i+2);
        %evaluate the ith distance constraint
        e_vec(i) = single_string_error_func(xA,yA,xB,yB,l_max);
    end
end