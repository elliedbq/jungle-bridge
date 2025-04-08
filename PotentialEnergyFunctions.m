function U_total = total_potential(coords, params)
    % Compute total potential energy of the bridge
    U_rb = PotentialEnergyFunctions.total_RB_potential(coords, params);
    U_g = PotentialEnergyFunctions.total_G_potential(coords, params);
    U_total = U_rb + U_g;
end

function U_RB_total = total_RB_potential(coords, params)
    % Compute total rubber band potential energy
    U_RB_total = 0;
    coords = [params.r0; coords; params.rn];
    
    for i = 1:params.num_links
        l0 = params.l0_list(i);
        k = params.k_list(i);
        xA = coords(2*i-1);
        yA = coords(2*i);
        xB = coords(2*i+1);
        yB = coords(2*i+2);
        
        U_RB_total = U_RB_total + PotentialEnergyFunctions.single_RB_potential(xA, yA, xB, yB, k, l0);
    end
end

function U_RB_i = single_RB_potential(xA, yA, xB, yB, k, l0)
    % Compute potential energy of a single rubber band
    l = sqrt((xB-xA)^2 + (yB-yA)^2);
    U_RB_i = (1/2)*k*(max(l - l0, 0))^2;
end

function U_G_total = total_G_potential(coords, params)
    % Compute total gravitational potential energy
    U_G_total = 0;
    
    for i = 1:params.num_links-1
        m = params.m_list(i);
        y = coords(2*i);
        U_G_total = U_G_total + PotentialEnergyFunctions.single_G_potential(y, m, params.g);
    end
end

function U_G_i = single_G_potential(y, m, g)
    % Compute gravitational potential energy of a single weight
    U_G_i = m * g * y;
end