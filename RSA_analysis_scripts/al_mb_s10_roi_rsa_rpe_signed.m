%% Run RSA for positive and negative PE trials in feedback ROIs

subInfo = readtable('Adult_Run_List.csv');
scriptPath = pwd; 

%% 2) Execute the roi searchlights

% IMPORTANT: specify rois here
roiPath = [pwd,'/ROI_combinedGLM','/feedback_identity']; 
roiList = dir(fullfile(roiPath,'*.nii'));

% IMPORTANT: specify which model you want to look at: 
% 1 = identity, 2 = qvals
RDM_roi = 1;
phase = 'feedback'; 

% IMPORTANT: specify your input directory here
output_folder = [pwd,'/CombinedGLM_Feedback_3mm_Multimodel_Searchlight_9mm_1ss'];

% define user options here 
userOptions = defineUserOptions_socialRisk(scriptPath);

for subIdx = 1:30

    % subjects that will be included in final analysis 
    if subInfo.exclusionCode(subIdx) == 0 
        
    % extract out the relevant onset info for both runs
    subNum = subInfo.subNum(subIdx); 
    
    % run for TG & SM
     for condition = 1:2
        % construct predictor RDMs for high RPE trials only 
        [pred_RDM,conditionVec] = construct_pred_RDMs_RSA_rpe_signed_feedback(subNum,condition,scriptPath,RDM_roi); % change this so you're only running relevant RDM
        extract_roi_coeffs_rpe(subNum,scriptPath,pred_RDM,userOptions,condition,conditionVec,roiPath,roiList,output_folder,phase)
     end 
    
    elseif subInfo.exclusionCode(subIdx) == 1
        fprintf('Invalid Subject'); 
    end 
    
end






