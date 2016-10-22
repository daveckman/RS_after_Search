function [true_means, obs_means, true_best_system_id, in_PZ] = RealSearchLog(k, common_var, n0)

% Maximize g(x) = ceil(log_2(x)) on the interval [1/16, 16];

% Std dev of sample mean based on n_0 samples.
samplemean_sd = sqrt(common_var/n0);

% Initialize k systems
true_means = zeros(1,k);
obs_means = zeros(1,k);

% Start first system at x = 0.75.
true_means(1) = ceil(log2(0.75));
obs_means(1) = normrnd(true_means(1), samplemean_sd);
best_system = 0.75;
best_system_id = 1;

for new_system_id = 2:k

    % Identify next system. Project onto [1/16, 16] if necessary.
    new_system = best_system + unifrnd(-1, 1);
    if new_system < 1./16;
        new_system = 1./16;
    elseif new_system > 16;
        new_system = 16;
    end

    % Generate sample mean based on n_0 samples of new system
    true_means(new_system_id) = ceil(log2(new_system));
    obs_means(new_system_id) = normrnd(true_means(new_system_id), samplemean_sd);
    
    % Update identity of best-looking system
    if obs_means(new_system_id) > obs_means(best_system_id);
        best_system = new_system;
        best_system_id = new_system_id;
    end
end

% Check if returned configuration is in the preference zone
in_PZ = 0;
[sorted_true_means, sorted_true_means_id] = sort(true_means);
if sorted_true_means(k) ~= sorted_true_means(k-1); % unique best system
    in_PZ = 1;
end
true_best_system_id = sorted_true_means_id(k);