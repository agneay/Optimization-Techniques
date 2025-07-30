clc; 
clear;
 
r4 = @(a) round(a, 4);

disp('Experiment-2. Implementation of Golden Section Search,')
disp('      Fibonacci Search for single variable optimization problems')
disp('=================================================================')

%% ----------------  USER INPUT  ----------------------------
func = input('Enter the function f(x): ');     % Example: @(x) (x-3)^2 + 1
L    = input('Enter the left  endpoint (L): ');
R    = input('Enter the right endpoint (R): ');
% Reuse L and R from Golden Section for Fibonacci
L_fib = L;
R_fib = R;
eps  = input('Enter the tolerance ε : ');
kmax = input('Enter the maximum number of stages: ');
fprintf('\n');
%% ========== GOLDEN SECTION SEARCH =========================
disp('One‑dimensional minimisation via the Golden‑Section Search method')
disp('=================================================================')

rho = round((3 - sqrt(5)) / 2,4);

x1 = r4(L + rho * (R - L));
x2 = r4(R - rho * (R - L));
f1 = r4(func(x1));
f2 = r4(func(x2));

hist  = zeros(kmax, 8);      % k | L | R | x1 | x2 | f1 | f2 | |R-L|
where = strings(kmax, 1);    % direction for k+1

k = 1;
hist(k,:) = [k, L, R, x1, x2, f1, f2, R - L];

while (hist(k,8) > eps) && (k < kmax)
    if f1 < f2
        R = x2;
        decision = "Left";    % L is fixed in next step
    else
        L = x1;
        decision = "Right";   % R is fixed in next step
    end

    len = R - L;
    x1 = r4(L + rho * (R - L));
    x2 = r4(R - rho * (R - L));
    f1 = r4(func(x1));
    f2 = r4(func(x2));
    where(k) = decision;

    k = k + 1;
    hist(k,:) = [k, L, R, x1, x2, f1, f2, len];
end

hist  = hist(1:k, :);
where = where(1:k);
where(end) = "";  % No direction for last row

xmin = r4((L + R) / 2);
fmin = r4(func(xmin));

% ----- Print Golden Section Table -----
fprintf('\nGolden Section Table\n');

fprintf('======================================================================================================================\n');
fprintf('%-3s │ %-10s │ %-10s │ %-10s │ %-10s │ %-10s │ %-10s │ %-10s │ %-15s|\n', ...
        'k','L','R','x1','x2','f(x1)','f(x2)','|R-L|','Fixed Dir');
fprintf('%-3s │ %-10s │ %-10s │ %-10s │ %-10s │ %-10s │ %-10s │ %-10s │ %-15s|\n', ...
        ' ',' ',' ',' ',' ',' ',' ',' ',' for k+1');
fprintf('======================================================================================================================\n');

rowfmt = '%-3d │ %-10.4f │ %-10.4f │ %-10.4f │ %-10.4f │ %-10.4f │ %-10.4f │ %-10.4f │ %-15s|\n';
for r = 1:k
    fprintf(rowfmt, hist(r,1:8), where(r));
end

fprintf('======================================================================================================================\n');

fprintf('\nThus, x_min ∈ [%.4f, %.4f];\n', L, R);
fprintf('Hence, x* = (%.4f + %.4f)/2 = %.4f\n', L, R, xmin);
fprintf('Function value  f(x*) = %.5f\n', fmin);

fprintf('\n');
fprintf('\n');
%% ========== FIBONACCI SEARCH (Stable and Matches Golden Table) ============
disp('One‑dimensional minimisation via the Fibonacci Search method')
disp('=============================================================')

inputStr = input('Enter the maximum number of stages in Fibonacci (or type "auto"): ', 's');

if strcmpi(inputStr, 'auto')
    % Auto-compute number of stages based on uncertainty range
    a0 = L_fib; 
    b0 = R_fib;
    final_range = eps;     % final uncertainty (you can assign your own value)
    ratio = final_range / (b0 - a0);
    % target reduction ratio

    % Start generating Fibonacci numbers
    F = [1, 1];  % don't call fibonacci_sequence yet
    N = 2;
    while (1)/F(N) > ratio
        F(N+1) = F(N) + F(N-1);
        N = N + 1;
    end
    kmax = N-1;   % Set number of stages
    fprintf('Auto-calculated number of stages: %d\n', kmax);
else
    % Convert string input to numeric
    kmax = str2double(inputStr);
    if isnan(kmax) || kmax <= 0
        error('Invalid number of stages. Enter a positive number or "auto".');
    end
end

% Now you can safely generate Fibonacci sequence
F = fibonacci_sequence(kmax + 2);  % this line comes after kmax is known


% ✅ Now that kmax is known, generate the Fibonacci sequence

% Generate enough Fibonacci numbers
F = fibonacci_sequence(kmax + 2);  % We need F_0 to F_{kmax+1}
n = kmax+1;

% Initial interior points

x1 = r4(L_fib + (F(n-2) / F(n)) * (R_fib - L_fib));
x2 = r4(L_fib + R_fib - x1);
f1 = r4(func(x1));
f2 = r4(func(x2));

rsl = {};  % Store iteration data

for k = 1:kmax

    len = r4(R_fib - L_fib);

    if k<kmax

        if f1 < f2
            fixed = "Left";
        else
            fixed= "Right";
        end

    else
        fixed=' ';
    end

    rsl(k,:) = {k, L_fib, R_fib, x1, x2, f1, f2, len,fixed};

    L_min=L_fib;
    R_min=R_fib;

    if f1 < f2
        R_fib = x2;
        x2 = x1;
        f2 = f1;
        x1 = r4(L_fib + R_fib - x2);
        f1 = r4(func(x1));
    else
        L_fib = x1;
        x1 = x2;
        f1 = f2;
        x2 = r4(L_fib + R_fib - x1);
        f2 = r4(func(x2));
    end
   end
%% ----- Print Table in textbook format (fraction ratio + clean columns) -----
fprintf('\nFibonacci Search Table\n');
fprintf('=================================================================================================================\n');
fprintf('| %-3s | %-10s | %-10s | %-10s | %-10s | %-10s | %-10s | %-10s | %-10s |\n', ...
        'k', 'F_{n−k}/', 'L', 'R', 'x1', 'x2', 'f(x1)', 'f(x2)', 'Fixed Dir');
fprintf('| %-3s | %-10s | %-10s | %-10s | %-10s | %-10s | %-10s | %-10s | %-10s |\n', ...
        ' ', ' F_{n−k+1}', ' ', ' ', ' ', ' ', ' ', ' ', 'for k+1');
fprintf('=================================================================================================================\n');

rowfmt = '| %-3d | %-10s | %-10.4f | %-10.4f | %-10.4f | %-10.4f | %-10.4f | %-10.4f | %-10s |\n';

for r = 1:size(rsl,1)
    k     = rsl{r,1};
    Lval  = rsl{r,2};
    Rval  = rsl{r,3};
    x1    = rsl{r,4};
    x2    = rsl{r,5};
    f1    = rsl{r,6};
    f2    = rsl{r,7};
    fixed = rsl{r,9};

    % Compute Fibonacci ratio string
    num_idx = n - k;
    den_idx = n - k + 1;
    if num_idx >= 1 && den_idx <= length(F)
        ratio_str = sprintf('%d/%d', F(num_idx), F(den_idx));
    else
        ratio_str = '—';
    end

    fprintf(rowfmt, k, ratio_str, Lval, Rval, x1, x2, f1, f2, fixed);
end

fprintf('=================================================================================================================\n');

% Final output summary
x_opt = r4((L_min + R_min) / 2);
f_opt = r4(func(x_opt));
fprintf('\nThus, x_min ∈ [%.4f, %.4f];\n', L_min, R_min);
fprintf('Hence, x* = (%.4f + %.4f)/2 = %.4f\n', L_min, R_min, x_opt);
fprintf('Function value f(x*) = %.5f\n', f_opt);

% === Fibonacci Generator Function ===
function F = fibonacci_sequence(n)
    F = zeros(1, n);
    F(1:2) = 1;
    for i = 3:n
        F(i) = F(i - 1) + F(i - 2);
    end
end
%%
%%   Example for reference
% % % % % Experiment-2. Implementation of Golden Section Search,
% % % % %       Fibonacci Search for single variable optimization problems
% % % % % =================================================================
% % % % % Enter the function f(x): @(x) x.^4 -14*x.^3 + 60*x.^2 - 70*x;
% % % % % Enter the left  endpoint (L): 0
% % % % % Enter the right endpoint (R): 2
% % % % % Enter the tolerance ε : 0.3
% % % % % Enter the maximum number of stages: 7
% % % % % 
% % % % % One‑dimensional minimisation via the Golden‑Section Search method
% % % % % =================================================================
% % % % % 
% % % % % Golden Section Table
% % % % % ───────────────────────────────────────────────────────────────────────────────────────────────────────────────────
% % % % % k   │ L          │ R          │ x1         │ x2         │ f(x1)      │ f(x2)      │ |R-L|      │ Fixed Dir         |
% % % % %     │            │            │            │            │            │            │            │  for k+1          |
% % % % % ───────────────────────────────────────────────────────────────────────────────────────────────────────────────────
% % % % % 1   │ 0.0000     │ 2.0000     │ 0.7640     │ 1.2360     │ -24.3608   │ -18.9596   │ 2.0000     │ Left              |
% % % % % 2   │ 0.0000     │ 1.2360     │ 0.4722     │ 0.7638     │ -21.0999   │ -24.3605   │ 1.2360     │ Right             |
% % % % % 3   │ 0.4722     │ 1.2360     │ 0.7640     │ 0.9442     │ -24.3608   │ -23.5931   │ 0.7638     │ Left              |
% % % % % 4   │ 0.4722     │ 0.9442     │ 0.6525     │ 0.7639     │ -23.8376   │ -24.3606   │ 0.4720     │ Right             |
% % % % % 5   │ 0.6525     │ 0.9442     │ 0.7639     │ 0.8328     │ -24.3606   │ -24.2879   │ 0.2917     │                   |
% % % % % ───────────────────────────────────────────────────────────────────────────────────────────────────────────────────
% % % % % 
% % % % % Thus, x_min ∈ [0.6525, 0.9442];
% % % % % Hence, x* = (0.6525 + 0.9442)/2 = 0.7984
% % % % % Function value  f(x*) = -24.36020
% % % % % One‑dimensional minimisation via the Fibonacci Search method
% % % % % =============================================================
% % % % % Enter the maximum number of stages in Fibonacci (or type "auto"): auto
% % % % % Auto-calculated number of stages: 5
% % % % % 
% % % % % Fibonacci Search Table
% % % % % =================================================================================================================
% % % % % | k   | F_{n−k}/   | L          | R          | x1         | x2         | f(x1)      | f(x2)      | Fixed Dir  |
% % % % % |     |  F_{n−k+1} |            |            |            |            |            |            | for k+1    |
% % % % % =================================================================================================================
% % % % % | 1   | 5/8        | 0.0000     | 2.0000     | 0.7500     | 1.2500     | -24.3398   | -18.6523   | Left       |
% % % % % | 2   | 3/5        | 0.0000     | 1.2500     | 0.5000     | 0.7500     | -21.6875   | -24.3398   | Right      |
% % % % % | 3   | 2/3        | 0.5000     | 1.2500     | 0.7500     | 1.0000     | -24.3398   | -23.0000   | Left       |
% % % % % | 4   | 1/2        | 0.5000     | 1.0000     | 0.7500     | 0.7500     | -24.3398   | -24.3398   | Right      |
% % % % % | 5   | 1/1        | 0.7500     | 1.0000     | 0.7500     | 1.0000     | -24.3398   | -23.0000   |            |
% % % % % =================================================================================================================
% % % % % 
% % % % % Thus, x_min ∈ [0.7500, 1.0000];
% % % % % Hence, x* = (0.7500 + 1.0000)/2 = 0.8750
% % % % % Function value f(x*) = -24.10520

