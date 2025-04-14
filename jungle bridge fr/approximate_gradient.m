% given a function and a point, approximates and returns the gradient for
% that function at that point
function G = approximate_gradient(fun, V)
    h = 1e-5;
    G = zeros(size(V));
    for i = 1:length(V)
        dV = zeros(size(V));
        dV(i) = h;
        G(i) = (fun(V + dV) - fun(V - dV)) / (2 * h);
    end
end
