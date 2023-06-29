%% Batch the whole brain RSA for feedback phase 

%read subject info table with counter-balance and exclusion data 
subInfo = readtable('Adult_Run_List.csv');
scriptPath = pwd; 

%% 1) Create individual ROI brain masks for each participant

% for subIdx = 1:30
%     % subjects that will be included in final analysis 
%     if subInfo.exclusionCode(subIdx) == 0 
%         
%     % create a brain map for each subject 
%     subNum = subInfo.subNum(subIdx); 
%     create_individual_masks(subNum,scriptPath); 
%     end 
%     
%  end 

%% 2) Execute the whole-brain searchlight 

% define user options here 
userOptions = defineUserOptions_socialRisk(scriptPath);

% define your step size (1 for whole brain, 30 for testing)
ss = 1; 

for subIdx = 1:30

    % subjects that will be included in final analysis 
    if subInfo.exclusionCode(subIdx) == 0 
        
    subNum = subInfo.subNum(subIdx); 
    
    % run for TG & SM
    for condition = 1:2

        % construct predictor RDMs
        [pred_RDM,conditionVec] = construct_pred_RDMs_feedback(subNum,condition,scriptPath);
        % create contrasts
        execute_searchlight_feedback(subNum,scriptPath,pred_RDM,userOptions,condition,ss,conditionVec)
    end 
    
    elseif subInfo.exclusionCode(subIdx) == 1
        fprintf('Invalid Subject'); 
    end 
    
end

