% functions for day 20
% simulates the string bridge rough constrained optimization. 

function string_sim()
    %load data
    fname = 'string_data.csv';
    fpath = './';
    my_table = readtable([fpath, fname]);
    data = table2array(my_table(1:7, 2:3));
    x_coord = data(:,1)' ;
    y_coord = 30 - data(:,2)';
    mass_weight = table2array(my_table(1:5, 6));
    string_length = table2array(my_table(1:6, 9));

    % initialize system parameters
        %param_struct.r0 = [x_0;y_0]: coordinates of leftmost vertex
        % param_struct.rn = [x_n;y_n]: coordinates of rightmost vertex
        % param_struct.num_links: number of links in bridge
        % param_struct.l0_list = [l_1;...;l_n]: list of segment lengths
        % param_struct.m_list = [m_1;...;m_(n-1)]: list of weight masses
        % param_struct.g = 9.8 m/sec^2: gravitational acceleration
    param = struct();
    param.r0 = [x_coord(1); y_coord(1)];
    param.rn = [x_coord(length(x_coord)); y_coord(length(y_coord))];
    param.num_links = length(string_length);
    param.l0_list = string_length;
    param.k_list = 10;
    param.m_list = mass_weight;
    param.g = 9.8;

    [x_cstr, y_cstr] = generate_shape_prediction_constrained(param);
    [x_uncstr, y_uncstr] = generate_shape_prediction_unconstrained(param);

    figure; hold on;
    plot(x_coord, y_coord, '.-')
    plot(x_cstr, y_cstr, '.--')
    plot(x_uncstr, y_uncstr, '.--')
    legend('Measured', 'Constrained', 'UNconstrained')


end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% helper functions


function e_len = single_string_error_func(xA,yA,xB,yB,l_max)
%computes the distance constraint error for a pair of vertices
    %INPUTS:
        %(xA, yA): coordinates of first vertex
        %(xB, yB): coordinates of second vertex
        %l_max: maximum allowable distance between two vertices
    %OUTPUTS:
        %e_len: constraint error.
            % e_len<=0 when constraint is satisfied
            % e_len>0 when constraint is violated
    e_len = sqrt((xB-xA)^2 + (yB-yA)^2) - l_max;
end

function [e_vec,dummy] = bridge_error_func(coords,param)
%evaluates the distance constraint error across all links
    %INPUTS:
        %coords: vector of vertex positions from i=1 to i=(n-1)
        % [x_1;y_1;...;x_(n-1),y_(n-1)]
        %param_struct: struct containing parameters of the bridge
    %OUTPUTS:
        %e_val = [e_len1; ... ; e_len_n]: the vector of distance constraint errors
        %dummy = []: empty vector used to satisfy fmincon syntax
    
    %initialize error vector
    e_vec = zeros(param.num_links,1);

    %initialize dummy output for fmincon equality constraints
    dummy = [];

    %add the first and last vertex positions to the coordinate list
    coords = [param.r0;coords;param.rn];

    %iterate through each rubber band link
    for i = 1:param.num_links
        %extract the ith segment length
        l_max = param.l0_list(i);

        %extract the coordinates of the string ends
        xA = coords(2*i-1);
        yA = coords(2*i);
        xB = coords(2*i+1);
        yB = coords(2*i+2);

        %evaluate the ith distance constraint
        e_vec(i) = single_string_error_func(xA,yA,xB,yB,l_max);
    end
end

function u_g_total = total_g_potential_func(coords, param)
% computes the total potential energy of the string bridge
% bc string is not stretchy, only the gravitational energy matters.
    % u_total = u_gravity = sum of u_g_i
    % u_g_i = gravity*y coord*mass
    u_g_total = 0;
    
    for i = 1:length(coords)/2
        m = param.m_list(i);
        g = param.g;
        y = coords(i*2);

        u_g_i = m*g*y;

        u_g_total = u_g_total + u_g_i;
    end
end

function U_RB_i = single_RB_potential_func(xA,yA,xB,yB,k,l0)
% find the potential energy of one rubber band (string section)
    %compute stretched length of rubber band
    l = sqrt((xB-xA)^2 + (yB-yA)^2);
    %compute potential energy (remember to use max function!)
    U_RB_i = (1/2)*k*(max(l - l0, 0))^2;
end

function U_RB_total = total_RB_potential_func(coords, var)
% find the total potential energy of the rubber band (string section)
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

function u_total = total_potential_func(coords, param)
% finds total potential energy, including the potential energy of the
% string. this is used for unconstrained optimization bc we dont have the
% constraints that make sure the string is the right length
    u_total = total_g_potential_func(coords, param);
    u_total = u_total + total_RB_potential_func(coords, param);
end

function coords_guess = guess_coords(param)
    % generate a guess for coordinate locations
    x0 = param.r0(1);
    y0 = param.r0(2);
    xn = param.rn(1);
    yn = param.rn(2);
    x_guess = linspace(x0,xn,param.num_links+1);
    y_guess = linspace(y0,yn,param.num_links+1);
    % x and y need to be combined into one vector so that we can have one
    % multi dimensional point to put into fmincon as the initial vector. 
    coords_guess = zeros(2*(param.num_links-1),1);
    for n = 1:(param.num_links-1)
        coords_guess(2*n-1,1) = x_guess(n+1);
        coords_guess(2*n,1) = y_guess(n+1);
    end
end

function [x_list, y_list] = generate_shape_prediction_constrained(param)
% use fmincon (constrained optimization) to generate a prediction for the shape of the bridge.
    %INPUTS
        %param
    %OUTPUTS
        % x_list: x coordinate guesses
        % y_list: y coordinate guesses
    
    % generate a guess for coordinate locations
    coords_guess = guess_coords(param);

    %define cost function
    f_cost = @(v_in) total_g_potential_func(v_in, param);

    %define constraint function
    f_cstr = @(v_in) bridge_error_func(v_in, param);

    %use fmincon to compute predicted locations. 
    coords_sol = fmincon(f_cost,coords_guess,[],[],[],[],[],[],f_cstr);

    %unpack result and combine with r0 and rn from param_struct
    %to generate list of positions, x_list and y_list
    V_list = [param.r0;coords_sol;param.rn];
    x_list = V_list(1:2:(end-1));
    y_list = V_list(2:2:end);
end

function [x_list, y_list] = generate_shape_prediction_unconstrained(param)
% uses gradient descent (unconstrained optimization) to predict the shape of the bridge
    %INPUTS
        %param
    %OUTPUTS
        % x_list: x coordinate guesses
        % y_list: y coordinate guesses  
% because we have no constraints, we use k and l0 to make sure the string
% doesnt stretch. in the case of the string, k is infinite (the cost 
% function heavily penalizes stretching) and l0 is the measured length. 
    %specify optimization parameters
    opt_params = struct();
    opt_params.beta = .5;
    opt_params.gamma = .9;
    opt_params.max_iter = 500;
    opt_params.min_gradient = 1e-7;
    
    % initialize string stiffness
    param.k_list = ones(param.num_links,1);

    % generate a guess for coordinate locations 
    coords_guess = guess_coords(param);

    %run gradient descent mutliple times while increasing
    %the spring stiffness after each iteration
    %the solution to the current iteration is the initial guess
    %for the previous iteration.
    for n = 1:10
        %use anonymous function syntax to define the cost func
        %define cost func as the total potential energy function
        %using the current values in param_struct
        f_cost = @(V_in) total_potential_func(V_in,param);

        %use gradient descent function to compute
        %the predicted vertex locations (for the current iteration)
        coords_guess = run_gradient_descent(f_cost,coords_guess,opt_params);

        %increase the stiffnesses by a factor of 10
        param.k_list = param.k_list*10;
    end

    %the solution is coordinates values at the last iteration
    coords_sol = coords_guess;

    %unpack result and combine with r0 and rn from param_struct
    %to generate list of positions, x_list and y_list
    V_list = [param.r0;coords_sol;param.rn];
    x_list = V_list(1:2:(end-1));
    y_list = V_list(2:2:end);
end