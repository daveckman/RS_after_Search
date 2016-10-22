function [true_means, obs_means, obs_vars, best_system] = AdvSearchRinott(k, delta, common_var, n0, h)

% Std dev of sample mean based on n_0 replications
samplemean_sd = sqrt(common_var/n0);
sd = sqrt(common_var);

% Initialize k systems
obs_means = zeros(1,k);
true_means = zeros(1,k);
obs_vars = zeros(1,k);

% Sample from first system
new_data = normrnd(0, sd, n0, 1);
obs_means(1) = mean(new_data);
obs_vars(1) = var(new_data);

% Initialize statistics of best system and best-looking system
max_mean = obs_means(1);
max_system = 1;
best_mean = 0;
best_system = 1;

% Start adding systems (to mimic a search)
for new_system = 2:k
    if max_system == best_system % The best "looks" best
        % Introduce system that is delta better than the true best
        true_means(new_system) = best_mean + delta;
        best_mean = true_means(new_system);
        best_system = new_system;
    else % A non-best system "looks" best
        % Introduce system that is delta worse than the true best
        true_means(new_system) = best_mean - delta;
    end

    % Take n_0 replications of new system and calculate statistics
    new_data = normrnd(true_means(new_system), sd, n0, 1);
    obs_means(new_system) = mean(new_data);
    obs_vars(new_system) = var(new_data);

    % Calculate overall number of required replications
    N = max(ceil(h^2*obs_vars(new_system)/delta^2),n0);
    N_new = N - n0;

    % If more replications are needed, take them and update sample mean
    if N_new > 0 
        new_obs_mean = normrnd(true_means(new_system), sqrt(common_var/N_new));
        obs_means(new_system) = (1/N)*(obs_means(new_system)*n0 + new_obs_mean*N_new);
    end

    % Does the new system "look" better than the previous best-looking system?
    if obs_means(new_system) > max_mean;
        max_mean = obs_means(new_system);
        max_system = new_system;
    end
end