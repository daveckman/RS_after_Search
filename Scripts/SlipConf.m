function [obs_means, obs_vars] = SlipConf(k, delta, common_var, n0);

samplemean_sd = sqrt(common_var/n0);

obs_means = zeros(1, k);
obs_vars = zeros(1, k);

obs_means(1) = normrnd(delta, samplemean_sd);
obs_vars(1) = common_var*chi2rnd(n0-1)/(n0-1);

for i = 2:k
	obs_means(i) = normrnd(0, samplemean_sd);
	obs_vars(i) = common_var*chi2rnd(n0-1)/(n0-1);
end