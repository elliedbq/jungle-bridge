function [k_all, l0_all] = RubberBandAnalysis()
    % Main function to estimate rubber band parameters
    [mass, length] = load_rubber_band_data();
    k_all = zeros(1,6);
    l0_all = zeros(1,6);
    
    for i = 1:6
        [A, Y] = prepare_rubber_band_data(mass, length, i);
        [k, l0] = fit_rubber_band_model(A, Y);
        k_all(i) = k;
        l0_all(i) = l0;
    end
end

function plot_line_of_best_fit(rb_num)
    % Plot line of best fit for a specific rubber band
    [mass, length] = load_rubber_band_data();
    [A, Y] = prepare_rubber_band_data(mass, length, rb_num);
    [~, ~, m, b] = fit_rubber_band_model(A, Y);
    
    figure;
    hold on;
    plot(Y, A(:,1), '.');
    x_val = 0:100:1600;
    y_val = (x_val - b)./m;
    plot(x_val, y_val);
    xlabel('Force (N)');
    ylabel('Rubber Band Length (m)');
    title(['Line of Best Fit for Rubber Band ' num2str(rb_num)]);
    legend('Data points', 'Line of best fit', 'Location', 'northwest');
    hold off;
end

function plot_surface_plot(rb_num)
    % Plot error surface for a specific rubber band
    [mass, length] = load_rubber_band_data();
    [A, Y] = prepare_rubber_band_data(mass, length, rb_num);
    [~, ~, m, b] = fit_rubber_band_model(A, Y);
    
    x_vals = A(:,1);
    y_vals = Y;
    m_range = linspace(m - 20000, m + 20000, 200);
    b_range = linspace(b - 5000, b + 5000, 200);
    
    [M, B] = meshgrid(m_range, b_range);
    Z = zeros(size(M));
    
    for i = 1:length(x_vals)
        xi = x_vals(i);
        yi = y_vals(i);
        Z = Z + (M .* xi + B - yi).^2;
    end
    
    figure;
    surf(M, B, Z);
    xlabel('Slope m');
    ylabel('Intercept b');
    zlabel('Cost E(m, b)');
    title(['Surface Plot for Rubber Band ' num2str(rb_num)]);
    colorbar;
    colormap parula;
    shading interp;
    view(45, 30);
end

function [mass, length] = load_rubber_band_data()
    % Load rubber band data from CSV
    fname = 'jungle_bridge_data.csv';
    my_table = readtable(fname);
    data = table2array(my_table(1:24, 2:3));
    mass = data(:,1)';
    length = data(:,2)' ./ 100; % convert to meters
end

function [A, Y] = prepare_rubber_band_data(mass, length, rb_num)
    % Prepare data for a specific rubber band
    g = 9.8;
    num = rb_num * 4;
    mass_vals = mass(num - 3 : num)' / 1000;
    length_vals = length(num - 3 : num)';
    
    force_vals = g .* mass_vals;
    A = [length_vals, ones(size(length_vals))];
    Y = force_vals;
end

function [k, l0, m, b] = fit_rubber_band_model(A, Y)
    % Fit linear model to rubber band data
    q = (A'*A)\(A'*Y);
    m = q(1);
    b = q(2);
    k = m;
    l0 = -b/m;
end