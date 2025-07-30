clc;
clear;
close all;
format short;
%% -------- Symbolic and user inputs --------
syms x;
disp('------------------------------------------------------------------------')
disp('                   Secant method                                        ')
disp('------------------------------------------------------------------------')
func_str = input('Enter the function f(x)  : ', 's');   % e.g., '0.65 - (0.75/(1+x^2)) - 0.65*x*atan(1/x)'
f        = str2func(['@(x)', func_str]);

lambda1  = input('Enter initial point: ');     % Step 1: λ₁ = A = 0
t0       = input('Enter initial step size t0 (e.g., 0.1): ');
ep       = input('Enter the tolerance ε   (e.g., 1e-2) : ');
nmax     = input('Enter maximum number of iterations   : ');

% Optional bracket inputs
A_str = input('Enter lower bound A (press Enter to auto) : ', 's');
B_str = input('Enter upper bound B (press Enter to auto) : ', 's');

%% -------- Derivative --------
fsym   = str2sym(func_str);
df_sym = simplify(diff(fsym, x));
fprintf('\nSymbolic derivative f''(x):\n');
disp(df_sym);

%% -------- Bracketing (user or automatic) --------
A_user = ~isempty(A_str);
B_user = ~isempty(B_str);

if A_user && B_user
    A = str2double(A_str);
    B = str2double(B_str);
    fA = double(subs(df_sym, x, A));
    fB = double(subs(df_sym, x, B));

    if fA >= 0
        error('Invalid input: f''(A) must be negative. f''(%.2f) = %.4f', A, fA);
    elseif fB <= 0
        error('Invalid input: f''(B) must be positive. f''(%.2f) = %.4f', B, fB);
    else
        fprintf('User-supplied bracket accepted: A = %.4f (f''=%.6f), B = %.4f (f''=%.6f)\n', A, fA, B, fB);
    end
else
    %% -------- Automatic bracketing --------
    fprintf('\nAuto-bracketing in progress...\n');
    A = lambda1;
    fA = double(subs(df_sym, x, A));
    if fA >= 0
        error('f''(A=0) must be negative for algorithm to proceed.');
    end
    t = t0;
    max_expand = 50;
    for i = 1:max_expand
        B = lambda1+t;
        fB = double(subs(df_sym, x, B));
        if fB > 0
            break;
        end
        t = 2 * t;
        if fB <= 0
        A=B;
        fA=fB;
        end
    end
    fprintf('\nAuto-bracket found: A = %.4f (f''=%.6f), B = %.4f (f''=%.6f)\n', A, fA, B, fB);
end

%% -------- Secant-based iterations --------
converged = false;
rsl = [];  % iteration log

for k = 1:nmax
    % Step 5: Compute new λ
    lambda = A - fA * (B - A) / (fB - fA);
    fp_lambda = double(subs(df_sym, x, lambda));
    
    % Store results
    rsl(k,:) = [k, A, B, fA, fB, lambda,fp_lambda, abs(fp_lambda)];

    % Step 6: Convergence test
    if abs(fp_lambda) <= ep
        converged = true;
        break;
    end

    % Step 7–8: Update bracket
    if fp_lambda < 0
        A = lambda; fA = fp_lambda;
    else
        B = lambda; fB = fp_lambda;
    end
end
%% -------- Results --------
Vars = {'k','A','B','f''(A)','f''(B)','x_k+1','f''(x_k+1)','|f''(x_k+1)|'};
ResultTable = array2table(rsl, 'VariableNames', Vars)

xopt = lambda;
fmin = double(f(xopt));  
fprintf('\n Optimal x = %.6f\n', xopt);
fprintf('Optimal f(x) = %.6f\n', fmin); % Replace fmin with -fmin for maximization problem

if ~converged
    fprintf('Warning: method did not converge within %d iterations.\n', nmax);
end
%% %% Example 
% % Enter the function f(x)  : x.^3 - 12.2*x.^2 + 7.45*x - sin(x)
% % Enter initial point: 7.6
% % Enter initial step size t0 (e.g., 0.1): 0.2
% % Enter the tolerance ε   (e.g., 1e-2) : 0.01
% % Enter maximum number of iterations   : 100
% % Enter lower bound A (press Enter to auto) : 7.6
% % Enter upper bound B (press Enter to auto) : 14