clc;
clear;
close all;
format short;
%% Symbolic and user inputs
syms x;

% User inputs
disp('------------------------------------------------------------------------')
disp('                   Newton method                                        ')
disp('------------------------------------------------------------------------')
func_str = input('Enter the function f(x): ', 's');  % e.g., '2*sin(x) - x.^2/10'
f = str2func(['@(x)', func_str]);

x0 = input('Enter the initial guess x0: ');           % e.g., 2.5
ep = input('Enter the tolerance (e.g., 1e-5): ');      % e.g., 1e-5
nmax = input('Enter the maximum number of iterations: ');  % e.g., 20
%% Derivatives using symbolic math
fsym = str2sym(func_str);
df_sym = diff(fsym, x);     % First derivative (symbolic)
d2f_sym = diff(df_sym, x);  % Second derivative (symbolic)

% Simplify expressions
df_sym = simplify(expand(df_sym));
d2f_sym = simplify(expand(d2f_sym));
%% Display symbolic derivatives
fprintf('\nFirst derivative f''(x):\n');
disp(df_sym);

fprintf('Second derivative f''''(x):\n');
disp(d2f_sym);
%% Initialize Newton's method
xx(1) = x0;
converged = false;

% Newton's method loop
for k = 1:nmax
    xk = xx(k);
    df_val = double(subs(df_sym, x, xk));
    d2f_val = double(subs(d2f_sym, x, xk));

    if d2f_val == 0
        error('Second derivative is zero. Newton''s method fails.');
    end

    x_next = xk - df_val / d2f_val;
    err = abs(x_next - xk);

    % Store results with derivative values
    rsl(k,:) = [k-1, xk, df_val, d2f_val, x_next, err];

    xx(k+1) = x_next;

    if err < ep
        converged = true;
        break;
    end
end
%% Table display
Variables = {'k', 'x(k)', 'f''(x)', 'f''''(x)', 'x(k+1)', 'Error'};
Resl = array2table(rsl, 'VariableNames', Variables)
%% Final output
xopt = xx(end);
fmin = f(xopt);  % Function is already numeric

fprintf('\nOptimal value of x = %.4f\n', xopt);
fprintf('Optimal value of f(x) = %.4f\n', -fmin);% Replace -fmin with fmin for minimization problem

if ~converged
    fprintf('Warning: Method did not converge within %d iterations.\n', nmax);
end
%% Example 1
% % % Enter the function f(x): 0.65-[0.75/(1+x.^2)] - 0.65*x*atan(1/x)
% % % Enter the initial guess x0: 0.1
% % % Enter the tolerance (e.g., 1e-5): 0.01
% % % Enter the maximum number of iterations: 100
% % % % 
% % % % First derivative f'(x):
% % % % (0.65*x)/(x^2 + 1) - 0.65*atan(1/x) + (1.5*x)/(x^2 + 1)^2
% % % % 
% % % % Second derivative f''(x):
% % % % -(0.4*(8.0*x^2 - 7.0))/(x^2 + 1.0)^3
% % % % 
% % % % 
% % % % Resl =
% % % % 
% % % %   4×6 table
% % % % 
% % % %     k     x(k)         f'(x)       f''(x)    x(k+1)       Error   
% % % %     _    _______    ___________    ______    _______    __________
% % % % 
% % % %     0        0.1       -0.74483    2.6866    0.37724       0.27724
% % % %     1    0.37724       -0.13823     1.573    0.46512      0.087879
% % % %     2    0.46512      -0.017907    1.1713    0.48041      0.015289
% % % %     3    0.48041    -0.00050348    1.1057    0.48086    0.00045536
% % % % 
% % % % 
% % % % Optimal value of x = 0.4809
% % % % Optimal value of f(x) = 0.3100
