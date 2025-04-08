%save in file named 'test_gradient_descent.m'
%routine to test your gradient descent routine
% % function test_gradient_descent()
%specify optimization parameters
opt_params = struct();
opt_params.beta = .5;
opt_params.gamma = .9;
opt_params.max_iter = 500;
opt_params.min_gradient = 1e-7;
%number of tests to perform
num_tests = 1000;
%iterate num_tests times
%for each iteration, create a new function f(X)
%and a new test input X0 to test the differentiator at
for n = 1:num_tests
%generate a randomized input dimension
input_dim = randi([1,15]);
%create a struct containing the test function parameters
test_params = struct();
%generate a matrix, test_params.A, with
%dimension input_dim x input_dim
%with each element element being randomly selected from
%a normal distribution
test_params.A = randn([input_dim,input_dim]);
%make test_params.A a positive definite matrix
test_params.A = test_params.A'*test_params.A+.5*eye(input_dim);
%generate a column vector, test_params.B, of length input_dim
%with each element element being randomly selected from
%a normal distribution
test_params.B = randn([input_dim,1]);
%create a new test function
%this is essentially a random second-order (quadratic) function
%with input dimension input_dim and output dimension output_dim
test_fun = @(X_in) optimization_test_function(X_in,test_params);
%compute the analytical global minimum of the test_function
Xopt_analytical = -test_params.A\test_params.B;
%define X0 as a column vector, of length input_dim
%with each element element being randomly selected from
%a normal distribution + the analytical minimum location
X0 = 10*randn(input_dim,1)+Xopt_analytical;
%compute the minimum location using
%your optimizer
X_opt_numerical = run_gradient_descent(test_fun,X0,opt_params);
%compare numerical and analytical solutions
%pull out the largest differing value between the two
largest_error = max(abs(X_opt_numerical-Xopt_analytical));
%if the analytical and numerical solutions differ by too much
%print 'fail!'
if largest_error>1e-5
disp('fail!');
end
end
% % end
%computes a quadratic function on input X
function f_val = optimization_test_function(X,params)
f_val = .5*X'*params.A*X+params.B'*X;
end
