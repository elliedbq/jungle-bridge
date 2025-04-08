%this function runs gradient descent to find the local minimum
%of an input function given some initial guess
    %INPUTS:
        %fun: the function we want to optimize
        %V0: the initial guess for gradient descent
        %params: a struct defining the optimization parameters
        %params.beta: threshold for choosing alpha (step size scaling)
        %params.gamma: growth/decay multiplier for backtracking line-search
        %params.max_iter: maximum number of iterations for gradient descent
        %params.min_gradient: termination condition for gradient descent
    %OUTPUTS:
        %Vopt: The guess for the local minimum of V0

function Vopt = run_gradient_descent(fun,v0,params)
    %unpack params
    beta = params.beta;
    gamma = params.gamma;
    max_iter = params.max_iter;
    min_gradient = params.min_gradient;

    % initialize variables
    alpha = 1;  % initial step size
    v = v0;     % current guess
    n = 0;      % number of iterations
    f = fun(v);  % current function value
    g = approximate_gradient(fun, v); % current gradient

    % loop until it converges or it reaches max_iterations
    while n< max_iter && norm(g)> min_gradient
        NG2 = norm(g)^2;

        %shrink alpha til the function decreases enough 
            % based on equation in notes.
        while fun(v-alpha*g) > f -beta *alpha*NG2
            alpha = alpha * gamma; %shrinks alpha by gamma
        end

        v = v - alpha*g; % update position

        % recompute gradient and function value for the next round
        f = fun(v);
        g = approximate_gradient(fun, v);
        n = n + 1;
        alpha = 1; % need to reset alpha at the end of each round so that it doesnt stay tiny
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
