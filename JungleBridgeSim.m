function [k_list, l0_list, predicted_coords, measured_coords] = JungleBridgeSim()
    % Main function that runs jungle bridge simulation
    
    % Load and analyze rubber band data
    [k_list, l0_list] = RubberBandAnalysis();
    
    % Load bridge configuration data
    [x_coord, y_coord, mass_weight] = loadBridgeData();
    
    % Initialize system parameters
    params = struct();
    params.r0 = [x_coord(1); y_coord(1)];
    params.rn = [x_coord(end); y_coord(end)];
    params.num_links = length(k_list);
    params.k_list = k_list';
    params.l0_list = l0_list';
    params.m_list = mass_weight';
    params.g = 9.8;
    
    % Compute predicted bridge shape
    [x_list, y_list] = generate_shape_prediction(params);
    
    % Package outputs
    predicted_coords = struct('x', x_list, 'y', y_list);
    measured_coords = struct('x', x_coord, 'y', y_coord);
end

function [x_coord, y_coord, mass_weight] = loadBridgeData()
    % Load bridge configuration data from CSV
    fname = 'jungle_bridge_data.csv';
    my_table = readtable(fname);
    data = table2array(my_table(1:7, 8:9));
    x_coord = data(:,1)';
    y_coord = 30 - data(:,2)';
    mass_weight = table2array(my_table(1:5, 12));
end

function [x_list, y_list] = generate_shape_prediction(params)
    % Predict bridge shape using gradient descent
    opt_params = struct();
    opt_params.beta = 0.5;
    opt_params.gamma = 0.9;
    opt_params.max_iter = 500;
    opt_params.min_gradient = 1e-7;
    
    % Call total_potential directly (it must be in the same file or on path)
    f_cost = @(V_in) total_potential(V_in, params);
    
    % Generate initial guess
    x_guess = linspace(params.r0(1), params.rn(1), params.num_links+1);
    y_guess = linspace(params.r0(2), params.rn(2), params.num_links+1);
    
    coords_guess = zeros(2*(params.num_links-1),1);
    for n = 1:(params.num_links-1)
        coords_guess(2*n-1,1) = x_guess(n+1);
        coords_guess(2*n,1) = y_guess(n+1);
    end
    
    coords_sol = run_gradient_descent(f_cost, coords_guess, opt_params);
    
    V_list = [params.r0; coords_sol; params.rn];
    x_list = V_list(1:2:(end-1));
    y_list = V_list(2:2:end);
end

% Include all potential energy functions in this file
function U_total = total_potential(coords, params)
    U_rb = total_RB_potential(coords, params);
    U_g = total_G_potential(coords, params);
    U_total = U_rb + U_g;
end

function U_RB_total = total_RB_potential(coords, params)
    U_RB_total = 0;
    coords = [params.r0; coords; params.rn];
    
    for i = 1:params.num_links
        l0 = params.l0_list(i);
        k = params.k_list(i);
        xA = coords(2*i-1);
        yA = coords(2*i);
        xB = coords(2*i+1);
        yB = coords(2*i+2);
        
        U_RB_total = U_RB_total + single_RB_potential(xA, yA, xB, yB, k, l0);
    end
end

function U_RB_i = single_RB_potential(xA, yA, xB, yB, k, l0)
    l = sqrt((xB-xA)^2 + (yB-yA)^2);
    U_RB_i = (1/2)*k*(max(l - l0, 0))^2;
end

function U_G_total = total_G_potential(coords, params)
    U_G_total = 0;
    
    for i = 1:params.num_links-1
        m = params.m_list(i);
        y = coords(2*i);
        U_G_total = U_G_total + single_G_potential(y, m, params.g);
    end
end

function U_G_i = single_G_potential(y, m, g)
    U_G_i = m * g * y;
end