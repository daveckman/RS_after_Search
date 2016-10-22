% RealSearch.m
% Plots Figure 3(b)

alpha = 0.05; % Maximum allowable Pr(Incorrect Selection)
delta = 1; % Indifference zone (IZ) parameter
common_var = 1; % Assumed known, common variance
n0 = 10; % Number of initial replications per system

% Read Bechhofer h in from text file
hB_array = textread('Bech_h_list.txt','%f');

R = 100000; % Number of macroreplications of the procedures

% # of systems (a range of values for a plot of PCS vs k)
k_settings = [2:2:10, 20:10:100];% 120:20:200];%, 250:50:500, 600:100:1000];
num_settings = length(k_settings);

% Initialize vector for storing PCS entries
PCS_RealSearch = zeros(1, num_settings); 
PGS_RealSearch = zeros(1, num_settings);

% Modified Gupta subset-selection procedure

for l = 1:length(k_settings)
    k = k_settings(l); % Fix number of systems
    h = hB_array(l); % Bechhofer h
    
    % Calculate yardstick for subset-selection comparisons
    yardstick = max(h*sqrt(2*common_var/n0) - delta,0);
    
    % Initialize vector for storing indicators of correct (and good) surival events and preferenze zone events (0
    % or 1)
    CS_RealSearch = zeros(1, R); 
    GS_RealSearch = zeros(1, R);
    PZ_RealSearch = zeros(1, R); 

    for r = 1:R
        % Run one simulation of realistic search for k systems
        [true_means, obs_means, true_best_system_id, in_PZ] = RealSearchLog(k, common_var, n0);

        PZ_RealSearch(r) = in_PZ;
        
        % If the best system survives, mark as correct selection
        if obs_means(true_best_system_id) >= max(obs_means) - yardstick
            CS_RealSearch(r) = 1;
        end
        
        % If a good system survives, mark as good selection       
        max_true_means = true_means(true_best_system_id);
        if max(obs_means(true_means == max_true_means)) >= max(obs_means) - yardstick
            GS_RealSearch(r) = 1;
        end 
    end
        
    % Record empirical PCS
    PCS_RealSearch(l) = sum(CS_RealSearch(PZ_RealSearch == 1))/sum(PZ_RealSearch);
    PGS_RealSearch(l) = sum(GS_RealSearch)/R;
    
    fprintf('Tested Modified Gupta after realistic search for k = %d. # of PZ instances was %d.\n', k, sum(PZ_RealSearch))
end

% Plot PCS vs # of systems
figure
plot(k_settings, PCS_RealSearch, '-ok', 'LineWidth', 2, 'MarkerFaceColor', 'black');
hold on
plot(k_settings, PGS_RealSearch, '--sk', 'LineWidth', 2, 'MarkerFaceColor', 'black');
hold off
legend('PCS | \mu in PZ(1)','PGS')
xlabel('No. of Returned Systems (k)', 'FontSize', 14);
ylabel('Probability of Selection Event', 'FontSize', 14);
title('Modified Gupta for Realistic Search', 'FontSize', 14);
