% functions for day 17: jungle bridge part 1 - finding k and l0
% and plots the line of best fit and contour plot for the regression
% regression cost function
% rubber band order: beige, yellow, green, black, blue, red

function [k_all, l0_all] = rubber_band()
% find k_all and l0_all
    [mass, length_rb] = load`_csv();
    k_all = zeros(1,6);
    l0_all = zeros(1,6);
    for i = 1:6 %num rubber bands
        [A, Y] = find_a_y(mass, length_rb, i);
        [k, l0, ~, ~] = find_k_l0(A, Y);
        k_all(i) = k;
        l0_all(i) = l0;
    end
    
% plot line of best fit and contour plot for rubber band 1
    %line_of_best_fit()
    %contour_plot()
end

function contour_plot()
% plots the contour plot for rubber band number 1
    rb_num = 1;
    [mass, length_rb] = load_csv();
    [A, Y] = find_a_y(mass, length_rb, rb_num);
    [~, ~, m_opt, b_opt] = find_k_l0(A, Y);

    force_list = Y;
    length_list = A(:, rb_num);

    delta_m = 1;
    delta_b = .1;

    m_range = linspace(m_opt-delta_m/2,m_opt+delta_m/2,101);
    b_range = linspace(b_opt-delta_b/2,b_opt+delta_b/2,101);

    [m_grid,b_grid] = meshgrid(m_range,b_range);

    cost_grid = zeros(size(m_grid)); %what we end up plotting
    for i = 1:length(force_list)
        xi = length_list(i);
        yi = force_list(i);
        cost_grid = cost_grid + (m_grid*xi+b_grid-yi).^2;
    end

    %choose the level sets of the contour plot so that they are spaced
    %quadratically between the min and max value on the grid
    min_val = min(min(cost_grid));
    max_val = max(max(cost_grid));
    
    dval = sqrt(max_val-min_val);
    levels = min_val + linspace(0,dval,25).^2;

    %generate contour plot
    figure();
    contourf(m_grid,b_grid,cost_grid,levels(1:end-1),'color','w');
    hold on
    plot(m_opt,b_opt,'ro','markerfacecolor','r','markersize',5);
    xlabel('m (newtons/meter)');
    ylabel('b (newtons)');
    title('Regression Cost Function for Rubber Band 1');
end


function line_of_best_fit()
% plots the line of best fit for rubber band number 1
    rb_num = 1;
    [mass, length] = load_csv();
    [A, Y] = find_a_y(mass, length, rb_num);
    [~, ~, m, b] = find_k_l0(A, Y);
    
    x_min = A(1,rb_num) - .01;
    x_max = A(end,rb_num) + .01;

    figure; hold on;
    plot(A(:,rb_num), Y,'.')
    x_val =x_min:.01:x_max;
    y_val = m.*x_val + b;
    plot(x_val, y_val)
    ylabel('Force (N)')
    xlabel('Rubber Band Length (m)')
    title('Line of Best Fit for Rubber Band 1')
    legend('data points', 'line of best fit', Location='northwest')
    hold off;
end

function [mass, length] = load_csv()
% loads the csv of the measured weights and rubber band lengths.
    fname = 'rubber_band_data.csv';
    fpath = './';
    my_table = readtable([fpath, fname]);
    data = table2array(my_table(1:24, 2:3));
    mass = data(:,1)';
    length = data(:,2)';
    length = length ./ 100; % convert to meters
    mass = mass ./ 1000; % convert to kg
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
% also returns m and b from y = mx + b (the line of best fit)
    q = (A'*A)\(A'*Y);
    m = q(1);
    b = q(2);

    k = m;
    l0 = -b/m;
end

