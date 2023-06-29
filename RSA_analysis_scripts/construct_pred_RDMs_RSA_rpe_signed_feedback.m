 function [pred_RDM,conditionVec] = construct_pred_RDMs_RSA_rpe_signed_feedback(subNum,condition,scriptPath,RDM_roi)
% construct predictor RDMs

% initialize a structure to store predictor RDMs
pred_RDM =  struct;

% get subject data
data = readtable("adult_behavData_2ca2lr.csv");
sub_dat = data(data.subjNum == subNum & data.condition == condition,:); 

% set up directories 
glmpath = [scriptPath,'/Full_DM_Timings_Combined/'];

% get the spm file for this subject
timePath = [scriptPath,'/Full_DM_Timings_Combined/'];
sub = char(sprintf("%d_timings", subNum));
spm_file = fullfile(timePath, sub, 'SPM.mat');
load(spm_file);

%% create conditionVec for trial number nested within playerNum 

% % zscore rpes 
sub_dat.delta_z = zscore(sub_dat.delta(:));

% can change the threshold here to check if binning artifact
neg_rpes = find(sub_dat.delta_z(:) < 0);
pos_rpes = find(sub_dat.delta_z(:) >= 0); 

% separate RPEs into pos and neg trials, 1 = pos, -1 = neg (not negative vs. positive)
rpe_sign = zeros(height(sub_dat),1);
rpe_sign(neg_rpes) = -1;
rpe_sign(pos_rpes) = 1; 
sub_dat.rpe_sign = rpe_sign; 

% specific GLM -> hypothesis matrix mapping &
% create a structure to store pos and neg RPE trials
conditionVec = struct; 
conditionVec.negRPE = zeros(length(SPM.xX.name),3);
conditionVec.posRPE = zeros(length(SPM.xX.name),3);

% for sm find out where trial 1 starts after motion regressors
sm_t1_idx_str = 'Sn(2) low_feedback_1*bf(1)';
sm_t1_idx = find(contains(SPM.xX.name,sm_t1_idx_str));

% intialize separate counters for tg & sm and neg vs. pos rpes  
if condition == 1
    neg_cnt = 0; 
    pos_cnt = 0; 
else
    neg_cnt = sm_t1_idx - 1;
    pos_cnt = sm_t1_idx - 1;
end 

for playerNum = 1:4

    for trial = 1:15

        % find the raw trial associated with this stimulus presentation
        trial_idx = find(sub_dat.playerNum == playerNum & sub_dat.binned_trial_num == trial);
        rpe_sign = sub_dat.rpe_sign(trial_idx);
        
        if rpe_sign == 1 % pos rpe
            pos_cnt = pos_cnt + 1; 
        else      
            neg_cnt = neg_cnt + 1; 
        end 

        if playerNum == 1
            trial_str = sprintf('Sn(%d) low_feedback_%d*bf(1)',condition,trial);
        elseif playerNum == 2
            trial_str = sprintf('Sn(%d) high_feedback_%d*bf(1)',condition,trial); 
        elseif playerNum == 3
            trial_str = sprintf('Sn(%d) neutral_feedback_%d*bf(1)',condition,trial); 
        elseif playerNum == 4
            trial_str = sprintf('Sn(%d) random_feedback_%d*bf(1)',condition,trial);
        end 

        if rpe_sign == 1 % pos rpe
            conditionVec.posRPE(pos_cnt,1) = find(contains(SPM.xX.name,trial_str));  
            conditionVec.posRPE(pos_cnt,2) = playerNum;
            conditionVec.posRPE(pos_cnt,3) = trial;
        else 
            conditionVec.negRPE(neg_cnt,1) = find(contains(SPM.xX.name,trial_str));
            conditionVec.negRPE(neg_cnt,2) = playerNum;
            conditionVec.negRPE(neg_cnt,3) = trial;
        end 

    end 

end 

%% construct identity RDM 

if RDM_roi == 1
    
   neg_rpe_identity_matrix = ones(length(neg_rpes),length(neg_rpes));
   pos_rpe_identity_matrix = ones(length(pos_rpes),length(pos_rpes));
   
   neg_rpe_dims = zeros(4,1);
   pos_rpe_dims = zeros(4,1);
   
   neg_cnt = 0; 
   pos_cnt = 0; 
  
    for playerNum = 1:4
        
        neg_dim = find(conditionVec.negRPE(:,2) == playerNum);
        pos_dim = find(conditionVec.posRPE(:,2) == playerNum);
     
        
        if playerNum == 1
            neg_cnt = neg_cnt + length(neg_dim);
            neg_rpe_identity_matrix(1:neg_cnt,1:neg_cnt) = 0 ;
            pos_cnt = pos_cnt + length(pos_dim);
            pos_rpe_identity_matrix(1:pos_cnt,1:pos_cnt) = 0 ;
            
        else
            neg_rpe_identity_matrix(neg_cnt+1:neg_cnt+length(neg_dim),neg_cnt+1:neg_cnt+length(neg_dim)) = 0;
            neg_cnt = neg_cnt + length(neg_dim);
            pos_rpe_identity_matrix(pos_cnt+1:pos_cnt+length(pos_dim),pos_cnt+1:pos_cnt+length(pos_dim)) = 0;
            pos_cnt = pos_cnt + length(pos_dim);
            
        end      
    end 
    
    % save
    pred_RDM.neg_rpe_identity_RDM = neg_rpe_identity_matrix; 
    pred_RDM.pos_rpe_identity_RDM = pos_rpe_identity_matrix; 
%     imagesc(pos_rpe_identity_matrix)
%     colorbar()

else 
    %% construct Qval RDM

    neg_rpe_qval_vec = ones(length(neg_rpes),1);
    pos_rpe_qval_vec = ones(length(pos_rpes),1);
    neg_cnt = 0;
    pos_cnt = 0;

    neg_rpe_qval_vec = []; 
    pos_rpe_qval_vec = []; 
    
    for playerNum = 1:4
            
        neg_rpe_idx = find(sub_dat.playerNum == playerNum & sub_dat.rpe_sign == -1);
        neg_rpe_qvals = sub_dat.q_val(neg_rpe_idx);
        neg_rpe_trial_nums = sub_dat.binned_trial_num(neg_rpe_idx);
        
        pos_rpe_idx = find(sub_dat.playerNum == playerNum & sub_dat.rpe_sign == 1);
        pos_rpe_qvals = sub_dat.q_val(pos_rpe_idx);
        pos_rpe_trial_nums = sub_dat.binned_trial_num(pos_rpe_idx);
        
        % initialize a temporary vector for player type
        temp_neg_rpe_qval = [];
        temp_pos_rpe_qval = [];
        
        if ~isempty(neg_rpe_trial_nums)
        
            % store trial 1 qvals for neg rpe trials
            if neg_rpe_trial_nums(1) == 1
                temp_neg_rpe_qval = [neg_rpe_qvals(1),temp_neg_rpe_qval];
            else 
                % get the previous q-val regardless of rpe type
                effective_trial_1 = neg_rpe_trial_nums(1) - 1; 
                prev_qval_idx = find(sub_dat.playerNum == playerNum & sub_dat.binned_trial_num == effective_trial_1);
                temp_neg_rpe_qval = [temp_neg_rpe_qval,sub_dat.q_val(prev_qval_idx)];   
            end 
        end 

        if ~isempty(pos_rpe_trial_nums)
            % store trial 1 qvals for pos rpe trials
            if pos_rpe_trial_nums(1) == 1
                temp_pos_rpe_qval = [pos_rpe_qvals(1),temp_pos_rpe_qval];
            else 
                % get the previous q-val regardless of rpe type
                effective_trial_1 = pos_rpe_trial_nums(1) - 1; 
                prev_qval_idx = find(sub_dat.playerNum == playerNum & sub_dat.binned_trial_num == effective_trial_1);
                temp_pos_rpe_qval = [temp_pos_rpe_qval,sub_dat.q_val(prev_qval_idx)];   
            end 
            
        end 

        % shift up q-vals by 1 trial
        neg_rpe_q_vals_shifted = neg_rpe_qvals(1:length(neg_rpe_qvals)-1);
        temp_neg_rpe_qval = [temp_neg_rpe_qval,neg_rpe_q_vals_shifted'];
        neg_rpe_qval_vec = [neg_rpe_qval_vec,temp_neg_rpe_qval];

        % shift up q-vals by 1 trial
        pos_rpe_q_vals_shifted = pos_rpe_qvals(1:length(pos_rpe_qvals)-1);
        temp_pos_rpe_qval = [temp_pos_rpe_qval,pos_rpe_q_vals_shifted'];
        pos_rpe_qval_vec = [pos_rpe_qval_vec,temp_pos_rpe_qval];
  
    end 

    % z-score the r_vec
    neg_rpe_qval_vec = zscore(neg_rpe_qval_vec); 
    [neg_rpe_qval_matrix] = construct_rdm(neg_rpe_qval_vec);
    
    pos_rpe_qval_vec = zscore(pos_rpe_qval_vec); 
    [pos_rpe_qval_matrix] = construct_rdm(pos_rpe_qval_vec);

    % save 
    pred_RDM.neg_rpe_qval_RDM = neg_rpe_qval_matrix;
    pred_RDM.pos_rpe_qval_RDM = pos_rpe_qval_matrix;
    
%     imagesc(neg_rpe_qval_matrix)s
%     colorbar()
%     imagesc(pos_rpe_qval_matrix)
%     colorbar()

end 
%% scratch

% length(fieldnames(pred_RDM));

end

