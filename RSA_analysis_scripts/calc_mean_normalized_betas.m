function [mean_betas] = calc_mean_normalized_betas(betas,conditionVec)

% this function will construct the neural RDM within our searchlight and
% extract the lower triangle to correlate with our hypothesis RDMs

% first full out regressors of interest & reoreganize betas to aligned with
% hypothesis RDM
r_interest = conditionVec(find(conditionVec ~= 0));
betas_ordered = betas(r_interest,:);
mean_betas = mean(betas_ordered');

end
