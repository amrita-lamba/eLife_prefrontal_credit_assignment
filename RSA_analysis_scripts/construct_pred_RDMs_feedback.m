function [pred_RDM,conditionVec] = construct_pred_RDMs_feedback(subNum,condition,scriptPath)
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

% specific GLM -> hypothesis matrix mapping 
conditionVec = zeros(length(SPM.xX.name),2);

% for sm find out where trial 1 starts after motion regressors
sm_t1_idx_str = 'Sn(2) low_feedback_1*bf(1)';
sm_t1_idx = find(contains(SPM.xX.name,sm_t1_idx_str));

% intialize counter 
if condition == 1
    cnt = 0; 
else
    cnt = sm_t1_idx - 1;
end 


for playerNum = 1:4
    
    for trial = 1:15
        
        % determine if subject responded on trial 
        t_data = sub_dat(sub_dat.playerNum == playerNum & sub_dat.binned_trial_num == trial,:);
        cur_choice = t_data.endowment; 
        
        if cur_choice ~= -1 
        
        cnt = cnt + 1; 
        
        if playerNum == 1
            trial_str = sprintf('Sn(%d) low_feedback_%d*bf(1)',condition,trial);
        elseif playerNum == 2
            trial_str = sprintf('Sn(%d) high_feedback_%d*bf(1)',condition,trial); 
        elseif playerNum == 3
            trial_str = sprintf('Sn(%d) neutral_feedback_%d*bf(1)',condition,trial); 
        elseif playerNum == 4
            trial_str = sprintf('Sn(%d) random_feedback_%d*bf(1)',condition,trial);
        end 
        
        conditionVec(cnt) = find(contains(SPM.xX.name,trial_str));   
        conditionVec(cnt,2) = playerNum;
        
        end 
    end 
end 

%% construct identity RDM 

% get number of responsded trials 
t_resp = length(unique(conditionVec(:,1)))-1;
identity_matrix = ones(t_resp,t_resp);  
   
cnt = 0; 
for playerNum = 1:4

    p_dim = find(conditionVec(:,2) == playerNum);

    if playerNum == 1
        cnt = cnt + length(p_dim);
        identity_matrix(1:cnt,1:cnt) = 0 ;

    else
        identity_matrix(cnt+1:cnt+length(p_dim),cnt+1:cnt+length(p_dim)) = 0;
        cnt = cnt + length(p_dim);

    end      
end 

% identity_matrix = ones(60,60);
% identity_matrix(1:15,1:15) = 0;
% identity_matrix(16:30,16:30) = 0;
% identity_matrix(31:45,31:45) = 0;
% identity_matrix(46:60,46:60) = 0; 

% save
pred_RDM.identity_RDM = identity_matrix; 
% imagesc(identity_matrix)
% colorbar()

%% construct Qval RDM

qval_vec = zeros(t_resp,1); 
cnt = 0;

for playerNum = 1:4

    p_trl = find(sub_dat.playerNum == playerNum);
    p_qval_val = sub_dat.q_val(p_trl); 

    for trial = 1:15
        
        % determine if subject responded on trial 
        t_data = sub_dat(sub_dat.playerNum == playerNum & sub_dat.binned_trial_num == trial,:);
        cur_choice = t_data.endowment; 
        
        if cur_choice ~= -1 
        
            cnt = cnt + 1; 

            if trial > 1
            % update vector
            qval_vec(cnt) = p_qval_val(trial); % shift by one trial (link BOLD image with qval on PREVIOUS trial)
            else
            qval_vec(cnt) = p_qval_val(trial);    
            end
        
        end 
    end 
end 

% z-score the r_vec
qval_vec = zscore(qval_vec); 
[qval_matrix] = construct_rdm(qval_vec);

% save
pred_RDM.qval_RDM = qval_matrix; 
% imagesc(qval_matrix)
% colorbar()

%% construct trial RDM

trial_vec = zeros(t_resp,1);
cnt = 0; 

for playerNum = 1:4

    for trial = 1:15
        
        % determine if subject responded on trial 
        t_data = sub_dat(sub_dat.playerNum == playerNum & sub_dat.binned_trial_num == trial,:);
        cur_choice = t_data.endowment; 
        
        if cur_choice ~= -1 
        
            cnt = cnt + 1; 
            trial_vec(cnt) = trial; 
            
        end 
    end 
end 

[trial_matrix] = construct_rdm(trial_vec);
pred_RDM.trial_RDM = trial_matrix; 
% imagesc(trial_matrix)
% colorbar()

end

