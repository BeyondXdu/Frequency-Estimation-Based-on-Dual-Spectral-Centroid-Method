function [T_n, p_value] = henzeZirklerTest(data, n_simulations, alpha)
% HENZEZIRKLERTEST Performs the multivariate normality test proposed by Henze and Zirkler (1990).
%
% Input:
%   data (N x p matrix): N observations of p-dimensional samples.
%   n_simulations (integer): Number of Monte Carlo simulations (recommended >= 1000).
%   alpha (scalar): Significance level (e.g., 0.05).
%
% Output:
%   T_n (scalar): The observed Henze-Zirkler test statistic.
%   p_value (scalar): P-value obtained via Monte Carlo simulation.

if nargin < 3
    alpha = 0.05;
end
if nargin < 2
    n_simulations = 1000;
end

[N, p] = size(data);
confidence_level = 1 - alpha; % Confidence Level (1-alpha)

if N < 20 || N < p + 1
    error('Sample size N must be at least 20 and greater than dimension p+1.');
end

%% --- 1. Parameter Estimation and Pre-calculation ---
mu = mean(data)';
Sigma = cov(data);

try
    Sigma_inv = inv(Sigma);
catch
    disp('Warning: Covariance matrix is singular. Adding a small regularization term.');
    Sigma = Sigma + 1e-6 * eye(p);
    Sigma_inv = inv(Sigma);
end

% Calculate Square Root Inverse of Sigma for Standardization
[V, D] = eig(Sigma);
D_sqrt_inv = diag(1 ./ sqrt(diag(D)));
Sigma_sqrt_inv = V * D_sqrt_inv * V';

Centered_Data = data - mu';
Z = Centered_Data * Sigma_sqrt_inv'; 

%% --- 2. Calculate Observed Henze-Zirkler Statistic T_n ---
T_n = calculate_HZ_statistic(Z, N, p);

%% --- 3. Monte Carlo Simulation for P-value ---
T_simulated = zeros(n_simulations, 1);
Z_sim_pool = randn(N * n_simulations, p);

parfor k = 1:n_simulations
    Z_sim = Z_sim_pool((k-1)*N + 1 : k*N, :);
    T_simulated(k) = calculate_HZ_statistic(Z_sim, N, p);
end

p_value = sum(T_simulated >= T_n) / n_simulations;

%% --- 4. Command Line Output (Accurate Statistical Conclusion) ---
disp(' ');
disp('=====================================================');
disp('  Henze-Zirkler Multivariate Normality Test Results  ');
disp('=====================================================');
fprintf('Sample Size (N): %d, Dimension (p): %d\n', N, p);
fprintf('Significance Level (alpha): %.2f, Confidence Level (1-alpha): %.0f%%\n', alpha, confidence_level * 100);
fprintf('Observed HZ Statistic (T_n): %.4f\n', T_n);
fprintf('Monte Carlo P-value: %.4f\n', p_value);

fprintf('\n--- Hypothesis Test Conclusion ---\n');
if p_value < alpha
    fprintf('Decision: Reject H0 (p < alpha). \nConclusion: The data is NOT multivariate normal.\n');
else
    % The "confidence" here refers to the level at which we accept the non-rejection.
    fprintf('Decision: Do NOT Reject H0 (p >= alpha). \nConclusion: The data CONFORMS to multivariate normal distribution at the %.0f%% confidence level.\n', confidence_level * 100);
end
disp('=====================================================');

%% --- 5. Visualization (Aesthetics Focus) ---
figure('Color', [1 1 1], 'Units', 'inches', 'Position', [1 1 6 5]);

% Histogram of the Simulated Distribution
h = histogram(T_simulated, 50, 'Normalization', 'pdf', ...
    'EdgeColor', 'none', 'FaceColor', [0.1 0.4 0.7], 'FaceAlpha', 0.8); 

hold on;

% Critical Value (Boundary of the Rejection Region)
T_crit = prctile(T_simulated, confidence_level * 100); 
plot([T_crit, T_crit], ylim, 'k--', 'LineWidth', 1.5, ...
    'DisplayName', [num2str(confidence_level * 100, '%.0f'), '%% Confidence Boundary']); 

% Observed T_n Value
plot([T_n, T_n], ylim, 'r-', 'LineWidth', 3, ...
    'DisplayName', ['Observed T_n = ', num2str(T_n, 4)]);

% Aesthetics
title('Monte Carlo Distribution of HZ Statistic (T_n)', 'FontSize', 14, 'FontWeight', 'bold');
xlabel('Test Statistic Value (T_n)', 'FontSize', 12);
ylabel('Probability Density', 'FontSize', 12);
legend('Location', 'northeast', 'FontSize', 10);
grid on;
box on; 
set(gca, 'TickDir', 'out', 'FontSize', 11);

hold off;

end

% --- Auxiliary Function: Calculate HZ Statistic T_n (Unchanged) ---
function T_n = calculate_HZ_statistic(Z, N, p)
% Calculates the HZ statistic based on standardized residuals Z (N x p).
% Z is assumed to be standardized, i.e., cov(Z) approx I.

a_sq = 0.5; % Smoothing parameter a^2 = 1/2

% 1. Calculate Squared Euclidean Distances D^2_j and D^2_jk
D_sq_j = sum(Z.^2, 2); % D^2_j (N x 1)
D_sq_jk_term1 = repmat(D_sq_j, 1, N);
D_sq_jk = D_sq_jk_term1 + D_sq_jk_term1' - 2 * (Z * Z'); % D^2_jk (N x N)

% 2. Term 1: Double Summation
Term1 = sum(sum(exp(-0.5 * a_sq * D_sq_jk)));

% 3. Term 2: Single Summation
Term2_base = sum(exp(-0.5 * a_sq / (1 + a_sq) * D_sq_j));
Term2 = 2 * (1 + a_sq)^(-p/2) * Term2_base;

% 4. Term 3: Constant Term
Term3 = N * (1 + 2 * a_sq)^(-p/2);

T_n = Term1 - Term2 + Term3;

end