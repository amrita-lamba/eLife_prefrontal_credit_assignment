function [pred_RDM,conditionVec] = construct_pred_RDMs_RSA_trialwise_choice(subNum,condition,scriptPath,RDM_roi)
% construct predictor RDMs

% initialize a structure to store predictor RDMs
pred_RDM =  struct;

% get subject data
data = readtable("adult_behavData_2ca2lr.xlsx");
sub_dat = data(data.subjNum == subNum & data.condition == condition,:); 

% get the spm file for this subject
timePath = [scriptPath,'/Full_DM_Timings_Combined/'];
sub = char(sprintf("%d_timings", subNum));
spm_file = fullfile(timePath, sub, 'SPM.mat');
load(spm_file);

%% create conditionVec for each trial number nested in the structure

% specific GLM -> hypothesis matrix mapping &
% create a structure to store high and low RPE trials
conditionVec = struct; 
for trial = 1:15

    trial_stuct = sprintf('trial_%d',trial);
    conditionVec.(trial_stuct) = zeros(length(SPM.xX.name),3);

end 

% for sm find out where trial 1 starts after motion regressors
sm_t1_idx_str = 'Sn(2) low_choice_1*bf(1)';
sm_t1_idx = find(contains(SPM.xX.name,sm_t1_idx_str));

for trial = 1:15
    
    % intialize separate counters for each trial
    if condition == 1
        cnt = 0; 
    else
        cnt = sm_t1_idx - 1;
    end
    
    for playerNum = 1:4

        % find the raw trial associated with this stimulus presentation
        trial_idx = find(sub_dat.playerNum == playerNum & sub_dat.binned_trial_num == trial);

        if playerNum == 1
            trial_str = sprintf('Sn(%d) low_choice_%d*bf(1)',condition,trial);
        elseif playerNum == 2
            trial_str = sprintf('Sn(%d) high_choice_%d*bf(1)',condition,trial); 
        elseif playerNum == 3
            trial_str = sprintf('Sn(%d) neutral_choice_%d*bf(1)',condition,trial); 
        elseif playerNum == 4
            trial_str = sprintf('Sn(%d) random_choice_%d*bf(1)',condition,trial);
        end
        
        % increment counter
        cnt = cnt + 1;
        trial_stuct = sprintf('trial_%d',trial);
        conditionVec.(trial_stuct)(cnt,1) = find(contains(SPM.xX.name,trial_str));  
        conditionVec.(trial_stuct)(cnt,2) = playerNum;
        conditionVec.(trial_stuct)(cnt,3) = trial_idx;

    end 
end 

%% construct identity RDM 

if RDM_roi == 1
    
    identity_matrix = [0,1,1,1,...
                      1,0,1,1,...
                      1,1,0,1,...
                      1,1,1,0];
    
    identity_matrix = reshape(identity_matrix,[4,4]);
    
    for trial = 1:15

        trial_stuct = sprintf('trial_%d',trial);
        pred_RDM.(trial_stuct) = identity_matrix;

    end 

else 
    %% construct Qval RDM
    
    for trial = 1:15
        
        qval_vec = zeros(4,1);
        
        for playerNum = 1:4
            
            % FIX THIS PART SO THAT THE PRED RDM INDEXES Q-VALUES FROM
            % PREVIOUS TRIAL 
            % find the raw trial associated with this stimulus presentation
            trial_idx = find(sub_dat.playerNum == playerNum & sub_dat.binned_trial_num == trial);
            qval_vec(playerNum) = sub_dat.q_val(trial_idx);     
            trial_stuct = sprintf('trial_%d',trial);
            pred_RDM.(trial_stuct) = identity_matrix;
        
        end 

    end 

end 

end

