function [corr_d] = calc_neural_corr_dist(betas,conditionVec)

% this function will construct the neural RDM within our searchlight and
% extract the lower triangle to correlate with our hypothesis RDMs

% first full out regressors of interest & reoreganize betas to aligned with
% hypothesis RDM
r_interest = conditionVec(find(conditionVec ~= 0));
betas_ordered = betas(r_interest,:);
dim = size(betas_ordered);
corr_distance_matrix = zeros(dim(1),dim(1));

% calculate correlation distance
for i = 1:dim(1) % each row 

    for j = 1:dim(1) % each column 

        % calculate correlation distance
        i_vec = betas_ordered(i,:);
        j_vec = betas_ordered(j,:);
        
        R = corrcoef(i_vec,j_vec);
        
        % populate matrix
        corr_distance_matrix(i,j) = 1 - R(2);
    
    end     

end

% imagesc(corr_distance_matrix)
% colorbar()

% extract the lower triangle
corr_d = corr_distance_matrix(find(tril(ones(size(corr_distance_matrix)),-1))); 

end
