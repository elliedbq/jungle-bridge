% functions for day 17: jungle bridge part 1 - finding k and l0
% rubber band order: beige, yellow, green, black, blue, red

function [k_all, l0_all] = day17()
% finds k and l0 for all rubber bands and returns them.
    [mass, length] = load_csv();
    k_all = zeros(1,6);
    l0_all = zeros(1,6);
    for i = 1:6 %num rubber bands
        [A, Y] = find_a_y(mass, length, i);
        [k, l0, ~, ~] = find_k_l0(A, Y);
        k_all(i) = k;
        l0_all(i) = l0;
    end

    line_of_best_fit(1)
    % contour_plot(mass, length, 1)
    surface_plot(mass, length, 1)

end

function line_of_best_fit(rb_num)
    [mass, length] = load_csv();
    [A, Y] = find_a_y(mass, length, rb_num);
    [~, ~, m, b] = find_k_l0(A, Y);

    figure; hold on;
    plot(Y,A(:,1), '.')
    x_val =0:100:1600;
    y_val = (x_val - b)./m;
    plot(x_val, y_val)
    xlabel('Force (N)')
    ylabel('Rubber Band Length (m)')
    title('line of best fit for rubber band 1')
    legend('data points', 'line of best fit', Location='northwest')
    hold off;
end

% function contour_plot(mass, length, rb_num)
%     [A, Y] = find_a_y(mass, length, rb_num);
% 
%     x_vals = A(:,1);
%     y_vals = Y;
% 
%     [m_vals, b_vals] = meshgrid(-50:1:50, -5:.01:5);
%     e_vals = zeros(size(m_vals));
% 
%     for i = 1:size(m_vals, 1)
%         for j = 1:size(m_vals, 2)
%             m = m_vals(i, j);
%             b = b_vals(i, j);
%             Y_pred = m * x_vals + b;
%             e_vals(i, j) = sum((Y_pred - y_vals).^2);
%         end
%     end
% 
%     figure;
%     contour(m_vals, b_vals, e_vals, 50)  % 50 contour lines
%     xlabel('Slope m')
%     ylabel('Intercept b')
%     title('Contour Plot of E(m, b)')
%     colorbar
% end

function surface_plot(mass, length, rb_num)
    [A, Y] = find_a_y(mass, length, rb_num);
    % length and force

    x_vals = A(:,1); % force
    y_vals = Y; % length
        [~, ~, m, b] = find_k_l0(A, Y);

    m_center = m;
    b_center = b;
    
    m_range = linspace(m_center - 20000, m_center + 20000, 200);  % 200 points
    b_range = linspace(b_center - 5000, b_center + 5000, 200);  

    [M, B] = meshgrid(m_range, b_range);
    Z = zeros(size(M));

    range = size(x_vals);
    disp(x_vals)
    % sum the errors for each layer of the mesh grid
    for i = 1:range
        xi = x_vals(i);
        yi = y_vals(i);
        Z = Z + (M .* xi + B - yi).^2;
    end

    figure();
    surf(M, B, Z);
    xlabel('Slope m');
    ylabel('Intercept b');
    zlabel('Cost E(m, b)')
    title(['Surface Plot of E(m, b) for Rubber Band ', num2str(rb_num)])
    colorbar
    colormap parula
    shading interp
    view(45, 30)

end

function [mass, length] = load_csv()
% loads the csv of the measured weights and rubber band lengths.
    fname = 'jungle_bridge_data.csv';
    fpath = './';
    my_table = readtable([fpath, fname]);
    data = table2array(my_table(1:24, 2:3));
    mass = data(:,1)';
    length = data(:,2)';
    length = length ./ 100; % convert to meters
end

function [A, Y] = find_a_y(mass, length, rubber_band_num)
% finds matrix A and Y for a single rubber band
    g = 9.8; %m/s^2
    num = rubber_band_num * 4;
    mass_vals = mass(num - 3 : num)';
    length_vals = length(num - 3 : num)';

    force_vals = g.*mass_vals; %compute force exerted by weight in Newtons
    A = [length_vals, ones(size(length_vals))];
    Y = force_vals;
end

function [k, l0, m, b] = find_k_l0(A, Y)
% finds k and l0 for individual rubber bands given A and Y
    q = (A'*A)\(A'*Y);
    m = q(1);
    b = q(2);

    k = m;
    l0 = -b/m;
end

