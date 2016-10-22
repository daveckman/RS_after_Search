function [true_means, obs_means, obs_vars, best_system] = AdvSearch(k, delta, common_var, n0, WantVar)

%delta = 1; % IZ parameter
%sigma = 1; % Standard deviation
%n0 = 10; % Initial sample sizes

samplemean_sd = sqrt(common_var/n0);
sd = sqrt(common_var);

% Initialize k systems
obs_means = zeros(1,k);
true_means = zeros(1,k);
obs_vars = zeros(1,k);

% Sample from first system
if strcmp(WantVar, 'N') == 1
    obs_means(1) = normrnd(0, samplemean_sd);
    %obs_vars(1) = common_var*chi2rnd(n0-1)/(n0-1);
else
    new_data = normrnd(0, sd, n0, 1);
    obs_means(1) = mean(new_data);
    obs_vars(1) = var(new_data);
end

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

    if strcmp(WantVar, 'N') == 1
        obs_means(new_system) = normrnd(true_means(new_system), samplemean_sd);
        %obs_vars(new_system) = common_var*chi2rnd(n0-1)/(n0-1);
    else
        new_data = normrnd(true_means(new_system), sd, n0, 1);
        obs_means(new_system) = mean(new_data);
        obs_vars(new_system) = var(new_data);
    end

    % Does the new system "look" better than the previous best-looking system?
    if obs_means(new_system) > max_mean;
        max_mean = obs_means(new_system);
        max_system = new_system;
    end
end