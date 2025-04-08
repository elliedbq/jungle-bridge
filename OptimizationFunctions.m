function [Vopt, V_history] = run_gradient_descent(fun, v0, params)
    % Run gradient descent optimization
    beta = params.beta;
    gamma = params.gamma;
    max_iter = params.max_iter;
    min_gradient = params.min_gradient;
    
    alpha = 1;
    v = v0;
    n = 0;
    f = fun(v);
    g = approximate_gradient(fun, v);
    
    V_history = v;
    
    while n < max_iter && norm(g) > min_gradient
        NG2 = norm(g)^2;
        
        while fun(v - alpha * g) > f - beta * alpha * NG2
            alpha = alpha * gamma;
        end
        
        v = v - alpha * g;
        V_history(:, end+1) = v;
        
        f = fun(v);
        g = approximate_gradient(fun, v);
        n = n + 1;
        alpha = 1;
    end
    
    Vopt = v;
end

function G = approximate_gradient(fun, V)
    h = 1e-5;
    G = zeros(size(V));
    for i = 1:length(V)
        dV = zeros(size(V));
        dV(i) = h;
        G(i) = (fun(V + dV) - fun(V - dV)) / (2 * h);
    end
end