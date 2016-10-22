% makePlots.m
% Generate Figures 1 and 2 in the paper.

alpha = 0.05; % Maximum allowable Pr(Incorrect Selection)
delta = 1; % Indifference zone (IZ) parameter
common_var = 1; % Assumed known, common variance
n0 = 10; % Number of initial replications per system

% Read Bechhofer h constant in from text file
hB_array = textread('Bech_h_list.txt','%f');
hR_array = textread('Rinott_h_list.txt','%f');

R = 10000; % Number of macroreplications of the procedures

% # of systems (a range of values for a plot of PCS vs k)
k_settings = [2:2:10, 20:10:100, 120:20:200, 250:50:500, 600:100:1000];
num_settings = length(k_settings);

% Initialize vectors for storing PCS entries

% Selection procedures

% Bechhofer
PCS_AS_Bech_n0 = zeros(1, num_settings);
PCS_AS_Bech_all = zeros(1, num_settings);
PCS_SC_Bech = zeros(1, num_settings);

% Rinott
PCS_AS_Rinott_n0 = zeros(1, num_settings);
PCS_AS_Rinott_all = zeros(1, num_settings);
PCS_SC_Rinott = zeros(1, num_settings);

% Subset-selection procedures

% Modified Gupta
PCS_AS_mod_Gupta = zeros(1, num_settings); 
PCS_SC_mod_Gupta = zeros(1, num_settings);

% Screen-to-the-Best
PCS_AS_STB = zeros(1, num_settings);
PCS_SC_STB = zeros(1, num_settings);

%%

% Bechhofer selection procedure

for l = 1:length(k_settings)
    k = k_settings(l); % Fix number of systems
    h = hB_array(l); % Bechhofer h
    Bech_n0 = ceil(2*h^2*common_var/delta^2);

    % Initialize vectors for correct selection events
    CS_AS_Bech_all = zeros(1, R); 
    CS_AS_Bech_n0 = zeros(1, R);
    CS_SC_Bech = zeros(1, R);

    for r = 1:R
        % Run one simulation of AS for k systems
        [true_means_AS, obs_means_AS, ~, best_system_AS] = AdvSearch(k, delta, common_var, Bech_n0, 'N');
        
        % If best system is chosen, mark as correct selection
        if max(obs_means_AS) == obs_means_AS(best_system_AS)
            CS_AS_Bech_all(r) = 1;
        end

        % Run one simulation of AS for k systems
        [true_means_AS, obs_means_AS, ~, best_system_AS] = AdvSearch(k, delta, common_var, n0, 'N');

        % Take additional replications of all systems
        N = max(ceil(2*h^2*common_var/delta^2),n0);
        N_new = N - n0;

        % If more replications are needed, update overall sample mean
        if N_new > 0 
            for i = 1:k
                new_obs_mean = normrnd(true_means_AS(i), sqrt(common_var/N_new));
                obs_means_AS(i) = (1/N)*(obs_means_AS(i)*n0 + new_obs_mean*N_new);
            end
        end

        % If best system is chosen, mark as correct selection
        if max(obs_means_AS) == obs_means_AS(best_system_AS)
            CS_AS_Bech_n0(r) = 1;
        end

        % Run one simulation in the slippage configuration
        [obs_means_SC, ~] = SlipConf(k, delta, common_var, Bech_n0);

        % If best system is chosen, mark as correct selection
        if max(obs_means_SC) == obs_means_SC(1)
            CS_SC_Bech(r) = 1;
        end
    end
    
    % Record empirical PCS
    PCS_AS_Bech_all(l) = sum(CS_AS_Bech_all)/R;
    PCS_AS_Bech_n0(l) = sum(CS_AS_Bech_n0)/R;
    PCS_SC_Bech(l) = sum(CS_SC_Bech)/R;

    fprintf('Tested Bechhofer after AS for k = %d.\n', k)
end

% Plot PCS vs # of systems
figure
plot(k_settings, PCS_AS_Bech_all, '-ok', 'LineWidth', 2, 'MarkerFaceColor', 'black');
hold on;
plot(k_settings, PCS_AS_Bech_n0, '-sk', 'LineWidth', 2, 'MarkerFaceColor', 'black');
plot(k_settings, PCS_SC_Bech, ':xk', 'LineWidth', 2);
hold off;

xlabel('No. of Returned Systems (k)', 'FontSize', 14);
ylabel('PCS', 'FontSize', 14);
title('Bechhofer', 'FontSize', 14);
legend('AS All', 'AS n_0', 'SC');

%%

% Rinott selection procedure

for l = 1:length(k_settings)
   k = k_settings(l); % Fix number of systems
   h = hR_array(l); % Bechhofer h
   
   % Initialize vector for storing correct surival events (1 or 0)
   CS_AS_Rinott_all = zeros(1,R);
   CS_AS_Rinott_n0 = zeros(1,R);
   CS_SC_Rinott = zeros(1,R);
   
   for r = 1:R
       % Run one simulation of AS for k systems
       [true_means_AS, obs_means_AS, obs_vars_AS, best_system_AS] = AdvSearchRinott(k, delta, common_var, n0, h);
       
       % If best system is chosen, mark as correct selection
       if max(obs_means_AS) == obs_means_AS(best_system_AS)
           CS_AS_Rinott_all(r) = 1;
       end

       % Run one simulation of AS for k systems
       [true_means_AS, obs_means_AS, obs_vars_AS, best_system_AS] = AdvSearch(k, delta, common_var, n0, 'Y');
       
       % Take additional replications of each systems
       for i = 1:k
           N = max(ceil(h^2*obs_vars_AS(i)/delta^2),n0);
           N_new = N - n0;

           % If more replications are needed, update overall sample mean
           if N_new > 0
               new_obs_mean = normrnd(true_means_AS(i), sqrt(common_var/N_new));
               obs_means_AS(i) = (1/N)*(obs_means_AS(i)*n0 + new_obs_mean*N_new);
           end
       end

       % If best system is chosen, mark as correct selection
       if max(obs_means_AS) == obs_means_AS(best_system_AS)
           CS_AS_Rinott_n0(r) = 1;
       end

       % Run one simulation in the slippage configuration
       [obs_means_SC, obs_vars_SC] = SlipConf(k, delta, common_var, n0);
       true_means_SC = zeros(1,k);
       true_means_SC(1) = delta;

       for i = 1:k
           N = max(ceil(h^2*obs_vars_SC(i)/delta^2),n0);
           N_new = N - n0;

           % If more replications are needed, update overall sample mean
           if N_new > 0
               new_obs_mean = normrnd(true_means_SC(i), sqrt(common_var/N_new));
               obs_means_SC(i) = (1/N)*(obs_means_SC(i)*n0 + new_obs_mean*N_new);
           end
       end

       % If best system is chosen, mark as correct selection
       if max(obs_means_SC) == obs_means_SC(1)
           CS_SC_Rinott(r) = 1;
       end
   end
   
   % Record empirical PCS
   PCS_AS_Rinott_all(l) = sum(CS_AS_Rinott_all)/R;
   PCS_AS_Rinott_n0(l) = sum(CS_AS_Rinott_n0)/R;
   PCS_SC_Rinott(l) = sum(CS_SC_Rinott)/R;

   fprintf('Tested Rinott after AS for k = %d.\n', k)
end

% Plot PCS vs # of systems
figure
plot(k_settings, PCS_AS_Rinott_all, '-ok', 'LineWidth', 2, 'MarkerFaceColor', 'black');
hold on;
plot(k_settings, PCS_AS_Rinott_n0, '-sk', 'LineWidth', 2, 'MarkerFaceColor', 'black');
plot(k_settings, PCS_SC_Rinott, ':xk', 'LineWidth', 2);
hold off;

xlabel('No. of Returned Systems (k)', 'FontSize', 14);
ylabel('Empirical PCS', 'FontSize', 14);
title('Rinott', 'FontSize', 14);
legend('AS All', 'AS n_0', 'SC');

%%

 % Modified Gupta subset-selection procedure
 
 for l = 1:length(k_settings)
     k = k_settings(l); % Fix number of systems
     h = hB_array(l); % Bechhofer h
     
     % Calculate yardstick for subset-selection comparisons
     yardstick = max(h*sqrt(2*common_var/n0) - delta,0);
     
     % Initialize vector for storing indicators of correct surival events (0 or 1)
     CS_AS = zeros(1, R); 
     CS_SC = zeros(1, R);
 
     for r = 1:R
         % Run one simulation of AS for k systems
         [true_means_AS, obs_means_AS, ~, best_system_AS] = AdvSearch(k, delta, common_var, n0, 'N');
 
         % If best system survives, mark as correct selection
         if obs_means_AS(best_system_AS) >= max(obs_means_AS) - yardstick
             CS_AS(r) = 1;
         end
 
         % Run one simulation in the slippage configuration
         [obs_means_SC, ~] = SlipConf(k, delta, common_var, n0);
 
         % If best system survives, mark as correct selection
         if obs_means_SC(1) >= max(obs_means_SC) - yardstick
             CS_SC(r) = 1;
         end
     end
     
     % Record empirical PCS
     PCS_AS_mod_Gupta(l) = sum(CS_AS)/R;
     PCS_SC_mod_Gupta(l) = sum(CS_SC)/R;
 
     fprintf('Tested Modified Gupta for k = %d.\n', k)
 end
 
 % Plot PCS vs # of systems
 figure
 plot(k_settings, PCS_AS_mod_Gupta, '-ok', 'LineWidth', 2, 'MarkerFaceColor', 'black');
 hold on;
 plot(k_settings, PCS_SC_mod_Gupta, ':xk', 'LineWidth', 2);
 hold off;
 xlabel('No. of Returned Systems (k)', 'FontSize', 14);
 ylabel('Empirical PCS', 'FontSize', 14);
 title('Modified Gupta', 'FontSize', 14);
 legend('AS','SC');

%%

% Screen-to-the-Best subset-selection procedure 

for l = 1:length(k_settings)
    k = k_settings(l); % Fix number of systems
    t = tinv((1-alpha)^(1/(k-1)),n0-1); % t = t_{(1-alpha)^(1/(k-1)),n0-1}
    display(t)
        
    % Initialize vector for storing indicators of correct surival events (0 or 1)
    CS_AS = zeros(1, R); 
    CS_SC = zeros(1, R);

    for r = 1:R
        % Run one simulation of AS for k systems
        [true_means_AS, obs_means_AS, obs_vars_AS, best_system_AS] = AdvSearch(k, delta, common_var, n0, 'Y');

        % Does the true best system survive subset selection for AS?
        W_bestj_AS = zeros(k, 1);
        for j = 1:k
            W_bestj_AS(j) = t*sqrt((obs_vars_AS(best_system_AS) + obs_vars_AS(j))/n0);
        end

        % If best systems survives, mark as correct selection
        if sum(obs_means_AS(best_system_AS) >= obs_means_AS' - max(W_bestj_AS - delta, 0)) == k;
            CS_AS(r) = 1;
        end
        
        % Run one simulation in the slippage configuration
        [obs_means_SC, obs_vars_SC] = SlipConf(k, delta, common_var, n0);

        % Do subset selection on system 1
        W_1j_SC = zeros(k,1);
        for j = 1:k;
            W_1j_SC(j) = t*sqrt((obs_vars_SC(1) + obs_vars_SC(j))/n0);
        end

        % If best systems survives, mark as correct selection
        if sum(obs_means_SC(1) >= obs_means_SC' - max(W_1j_SC - delta, 0)) == k;
            CS_SC(r) = 1;
        end
    end
    
    % Record empirical PCS
    PCS_AS_STB(l) = sum(CS_AS)/R;
    PCS_SC_STB(l) = sum(CS_SC)/R;

    fprintf('Tested Screen-to-the-Best for k = %d.\n', k)
end

% Plot PCS vs # of systems
figure
plot(k_settings, PCS_AS_STB, '-ok', 'LineWidth', 2, 'MarkerFaceColor', 'black');
hold on;
plot(k_settings, PCS_SC_STB, ':xk', 'LineWidth', 2);
hold off;
xlabel('No. of Returned Systems (k)', 'FontSize', 14);
ylabel('PCS', 'FontSize', 14);
title('Screen-to-the-Best', 'FontSize', 14);
legend('AS','SC');