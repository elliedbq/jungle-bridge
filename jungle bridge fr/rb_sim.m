% rubber band simulation

% main function that runs jungle bridge simulation code for rubber band:
% finds and plots the predicted rubber band bridge based on given
% parameters.
function rb_sim()
    % load data
    [k_list, l0_list] = rubber_band();
    fname = 'rubber_band_data.csv';
    fpath = './';
    my_table = readtable([fpath, fname]);
    data = table2array(my_table(1:7, 8:9));
    x_coord = data(:,1)';
    y_coord = 30 - data(:,2)';
    mass_weight = table2array(my_table(1:5, 12));

    %initialize the system parameters
    %which contains parameters describing behavior/measurements of bridge
        % param_struct.r0 = [x_0;y_0]: coordinates of leftmost vertex
        % param_struct.rn = [x_n;y_n]: coordinates of rightmost vertex
        % param_struct.num_links: number of rubber bands in bridge
        % param_struct.k_list = [k_1;...;k_n]: list of stiffnesses
        % param_struct.l0_list = [l0_1;...;l0_n]: list of natural lengths
        % param_struct.m_list = [m_1;...;m_(n-1)]: list of weight masses
        % param_struct.g = 9.8 m/sec^2: gravitational acceleration
        last = length(k_list);
    var = struct();
    var.r0 = [x_coord(1); y_coord(1)];
    var.rn = [x_coord(last+1); y_coord(last+1)];
    var.num_links = last;
    var.k_list = k_list';
    var.l0_list = l0_list';
    var.m_list = mass_weight';
    var.g = 9.8;

    %%%%%plots%%%%%
    % generate a plot comparing the predicted and measured bridge shape
    plot_predicted(var, x_coord, y_coord)

    % plot the progression of the predicted bridge shape
    plot_progression(var)
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% helper functions

%compute the potential energy of a SINGLE rubber band
    %INPUTS
        %(xA,yA): coordinates of left end of rubber band
        %(xB,yB): coordinates of right end of rubber band
        %k: stiffness of rubber band
        %l0: natural length of rubber band
    %OUTPUTS:
        %%U_RB_i: potential energy of rubber band
function U_RB_i = single_RB_potential_func(xA,yA,xB,yB,k,l0)
    %compute stretched length of rubber band
    l = sqrt((xB-xA)^2 + (yB-yA)^2);
    %compute potential energy (remember to use max function!)
    U_RB_i = (1/2)*k*(max(l - l0, 0))^2;
end

%compute the total potential energy of all rubber bands in bridge
    %INPUTS
        %param_struct: various variables
    %OUTPUTS
        %U_RB_total: total potential energy of all rubber bands
function U_RB_total = total_RB_potential_func(coords, var)
    U_RB_total = 0;
    coords = [var.r0; coords; var.rn];
    

    for i = 1:var.num_links
        l0 = var.l0_list(i);
        k = var.k_list(i);
        xA = coords(2*i-1); %left vertex of rubber band
        yA = coords(2*i);
        xB = coords(2*i+1); %right vertex of rubber band
        yB = coords(2*i+2);
        % used to find the length of a rubber band.

        U_RB_i = single_RB_potential_func(xA,yA,xB,yB,k,l0);
        
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
function U_g_i = single_G_potential_func(y, m, g)
    U_g_i = m*g*y;
end

%compute the total gravitational potential energy
%of all weights in bridge
%INPUTS:
    %var: various variables
%OUTPUTS
    %U_g_total: total gravitational energy
function U_g_total = total_G_potential_func(coords, var)
    U_g_total = 0;

    for i = 1:var.num_links-1
        m = var.m_list(i);
        y = coords(2*i);

        U_g_i = single_G_potential_func(y, m, var.g);
        
        U_g_total = U_g_total + U_g_i;
    end
end

%compute the total potential energy of the bridge
    %INPUTS
        %var: various variables
    %OUTPUTS
        %u_total: total potential energy of bridge
function u_total = total_potential_func(coords,var)
    u_g_total = total_G_potential_func(coords, var);
    u_rb_total = total_RB_potential_func(coords, var);
    u_total = u_rb_total + u_g_total;
end


%use gradient descent to predict the shape of the bridge
    % INPUTS
        % var: various variables
    %OUTPUTS
        % x_list = [x_0;x_1;...;x_n]: x coordinates of predicted vertex positions
        % y_list = [y_0;y_1;...;y_n]: y coordinates of predicted vertex positions
function [x_list, y_list] = generate_shape_prediction(param_struct, max_num)
    %specify optimization parameters
    opt_params = struct();
    opt_params.beta = .5;
    opt_params.gamma = .9;
    opt_params.max_iter = max_num;
    opt_params.min_gradient = 1e-7;

    [f_cost, coords_guess] = gradient_info(param_struct);

    coords_sol = run_gradient_descent(f_cost, coords_guess, opt_params);
    % coords_guess is a multi dimensional vector. so we can run gradient 
    % descent all at once.

    V_list = [param_struct.r0;coords_sol;param_struct.rn];
    x_list = V_list(1:2:(end-1));
    y_list = V_list(2:2:end);
end

% returns the cost function (total potential energy) and the vector
% containing the guess of coordinates for the rubber band bridge.
function [f_cost, coords_guess] = gradient_info(param_struct)
    % define cost function (the function we have to find the optimization of)
    f_cost = @(V_in) total_potential_func(V_in,param_struct);

    %generate an initial guess for the coordinate locations
    x0 = param_struct.r0(1);
    y0 = param_struct.r0(2);
    xn = param_struct.rn(1);
    yn = param_struct.rn(2);
    x_guess = linspace(x0,xn,param_struct.num_links+1); % [x0, x1, x2, ..., xn]
    y_guess = linspace(y0,yn,param_struct.num_links+1);

    coords_guess = zeros(2*(param_struct.num_links-1),1);
    % [x1;y1;x2;y2;xn;yn]
    for n = 1:(param_struct.num_links-1)
        coords_guess(2*n-1,1) = x_guess(n+1);
        coords_guess(2*n,1) = y_guess(n+1);
    end
end

%getting the second figure with the different alpha values. 
function plot_progression(param_struct)
    % initialize cost function and the guessed coordinates for the rubber
    % band bridge.

    % what the max_iterations will be
    iter_list = [1, 10, 25, 50, 100];

    figure(); hold on;
    for i = iter_list
        [x_vals, y_vals] = generate_shape_prediction(param_struct, i);
        plot(x_vals, y_vals, '.-')
    end
    legend('i = 1', 'i = 10', 'i = 25', 'i = 50', 'i = 100', location = 'southeast')
    xlabel('x direction (cm)')
    ylabel('y direction (cm)')
    title('Progression of Predictions')
    hold off
end

function plot_predicted(var, x_coords, y_coords)
% plots the predicted bridge vs. the measured bridge. 

    %compute the predicted bridge shape, max_iter = 500
    [x_list,y_list] = generate_shape_prediction(var, 500);

    %generate a plot comparing the predicted and measured bridge shape
    figure; hold on;
    plot(x_list, y_list, '.--')
    plot(x_coords, y_coords, '.-')
    legend('predicted', 'measured', location = 'southeast')
    xlabel('x direction (cm)')
    ylabel('y direction (cm)')
    title('Comparison Plot of Jungle Bridge')

end