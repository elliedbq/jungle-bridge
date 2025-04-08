%Main function that runs the jungle bridge simulation code
function jungle_bridge_sim()
    %initialize the system parameters
    %which contains parameters describing behavior/measurements of bridge
    % param_struct.r0 = [x_0;y_0]: coordinates of leftmost vertex
    % param_struct.rn = [x_n;y_n]: coordinates of rightmost vertex
    % param_struct.num_links: number of rubber bands in bridge
    % param_struct.k_list = [k_1;...;k_n]: list of stiffnesses (N/m)
    % param_struct.l0_list = [l0_1;...;l0_n]: list of natural lengths(meters)
    % param_struct.m_list = [m_1;...;m_(n-1)]: list of weight masses(kg)
    % param_struct.g = 9.8 m/sec^2: gravitational acceleration(m/s^2)
    param_struct = struct();
    x = [0
6.5
13.5
16.4
20
25.5
29.5];
    y = 20 - [3
12.5
17
17.75
16
10
2.5];
    param_struct.r0 = [x(1:length(x)-1),y(1:length(x)-1)];
    param_struct.rn = [x(2:length(x)),y(2:length(x))];
    param_struct.num_links = 6;
    param_struct.k_list = [39.92;121.01;240.13;39.45;137.83;107.84];
    param_struct.l0_list = [0.082289;0.078509;0.020779;0.025774;0.073011;0.078754];
    param_struct.m_list = [0.026;0.052;0.051;0.053;0.027];
    param_struct.g = 9.8;
    %compute the predicted bridge shape
    [x_list,y_list] = generate_shape_prediction(param_struct);
    %generate a plot comparing the predicted and measured bridge shape
    figure; hold on;
    plot(x_list,y_list)
    plot(x, y)
    legend('predicted', 'measured')
end

%compute the potential energy of a SINGLE rubber band
    %INPUTS
        %(xA,yA): coordinates of left end of rubber band
        %(xB,yB): coordinates of right end of rubber band
        %k: stiffness of rubber band
        %l0: natural length of rubber band
    %OUTPUTS:
        %U_RB_i: potential energy of rubber band
function U_RB_i = single_RB_potential_func(xA,yA,xB,yB,k,l0)
    %compute stretched length of rubber band
    l = sqrt((xA - xB)^2 + (yA - yB)^2);
    %compute potential energy (remember to use max function!)
    % aka 1/2 k x^2
    U_RB_i = (1/2)*k*(max((l - l0), 0))^2;
end

%compute the total potential energy of all rubber bands in bridge
    %INPUTS:
        %coords: vector of vertex positions from i=1 to i=(n-1)
        % [x_1;y_1;...;x_(n-1);y_(n-1)]
        %param_struct: struct containing parameters of the bridge
    %OUTPUTS:
        %U_RB_total: total potential energy of rubber bands in bridge
function U_RB_total = total_RB_potential_func(coords,param_struct)
    %initialize total spring potential energy
    U_RB_total = 0;
    %add the first and last vertex positions to the coordinate list
    disp(param_struct.r0)
    disp(coords)
    disp(param_struct.rn)
    coords = [param_struct.r0;coords;param_struct.rn];
    %iterate through each rubber band link
    for i = 1:param_struct.num_links
        %extract the ith stiffness and natural length
        l0 = param_struct.l0_list(i);
        k = param_struct.k_list(i);
        %extract the coordinates of the rubber band ends
        xA = coords(2*i-1);
        yA = coords(2*i);
        xB = coords(2*i+1);
        yB = coords(2*i+2);

        %compute the potential energy of the ith rubber band
        U_RB_i = single_RB_potential_func(xA,yA,xB,yB,k,l0);
        %add the ith potential to the total
        U_RB_total = U_RB_total + U_RB_i;
    end
end

%compute the gravitational potential energy of a SINGLE weight
    %INPUTS
        %(x,y): coordinates of a the current vertex
        %m: mass of the weight
        %g: gravitational acceleration
    %OUTPUTS:
        %U_g_i: gravitational potential energy of weight
function U_g_i = single_G_potential_func(x,y,m,g)
    %compute gravitational potential energy of weight
    U_g_i = m*g*y;
end

%compute the total gravitational potential energy
%of all weights in bridge
    %INPUTS:
        %coords: vector of vertex positions from i=1 to i=(n-1)
        % [x_1;y_1;...;x_(n-1),y_(n-1)]
        %param_struct: struct containing parameters of the bridge
    %OUTPUTS:
        %U_g_total: total gravitational potential energy
function U_g_total = total_G_potential_func(coords,param_struct)
    %initialize total gravitational potential energy
    U_g_total = 0;
    %iterate through each weight
    for i = 1:(param_struct.num_links-1)
        %extract the coordinates of the ith vertex
        x = coords(2*i-1);
        y = coords(2*i);
        %extract the ith mass and the gravitational acceleration
        m = param_struct.m_list(i);
        g = param_struct.g;
        %compute the gravitational potential energy of the ith mass
        U_g_i = single_G_potential_func(x,y,m,g);
        %add the ith potential to the total
        U_g_total = U_g_total + U_g_i;
    end
end

%compute the total potential energy of the bridge
    %INPUTS:
        %coords: vector of vertex positions from i=1 to i=(n-1)
        % [x_1;y_1;...;x_(n-1),y_(n-1)]
        %param_struct: struct containing parameters of the bridge
    %OUTPUTS:
        %U_total: total potential energy of the bridge
function U_total = total_potential_func(coords,param_struct)
    %compute the gravitational potential energy of the weights
    U_g_total = total_G_potential_func(coords,param_struct);
    %compute the spring potential energy of the rubber bands
    U_RB_total = total_RB_potential_func(coords,param_struct);
    %sum the two results
    U_total = U_g_total + U_RB_total;
end

%use gradient descent to predict the shape of the bridge
    %INPUTS:
        %param_struct: struct containing parameters of the bridge
    %OUTPUTS:
        %x_list = [x_0;x_1;...;x_n]: x coordinates of predicted vertex positions
        %y_list = [y_0;y_1;...;y_n]: x coordinates of predicted vertex positions
function [x_list,y_list] = generate_shape_prediction(param_struct)
    %specify optimization parameters
    opt_params = struct();
    opt_params.beta = .5;
    opt_params.gamma = .9;
    opt_params.max_iter = 500;
    opt_params.min_gradient = 1e-7;
    %use anonymous function syntax to define the cost func
    %define cost func as the total potential energy function
    %using the current values in param_struct
    f_cost = @(V_in) total_potential_func(V_in,param_struct);
    disp(V_in)
    %generate an initial guess for the coordinate locations
    %coords_guess = [x_1;y_1;...;x_(n-1);y_(n-1)]
    %I chose to use evenly spaced points from (x_0;y_0) to (x_n;y_n)
    x0 = param_struct.r0(1);
    y0 = param_struct.r0(2);
    xn = param_struct.rn(1);
    yn = param_struct.rn(2);
    x_guess = linspace(x0,xn,param_struct.num_links+1);
    y_guess = linspace(y0,yn,param_struct.num_links+1);
    coords_guess = zeros(2*(param_struct.num_links-1),1);
    for n = 1:(param_struct.num_links-1)
        coords_guess(2*n-1,1) = x_guess(n+1);
        coords_guess(2*n,1) = y_guess(n+1);
    end
    %use gradient descent function to compute
    %the predicted vertex locations
    coords_sol = run_gradient_descent(f_cost,coords_guess,opt_params);
    %unpack result and combine with r0 and rn from param_struct
    %to generate list of positions, x_list and y_list
    V_list = [param_struct.r0;coords_sol;param_struct.rn];
    x_list = V_list(1:2:(end-1));
    y_list = V_list(2:2:end);
end