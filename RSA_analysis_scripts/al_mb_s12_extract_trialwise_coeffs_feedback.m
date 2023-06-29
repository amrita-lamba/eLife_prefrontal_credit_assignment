%% Extract trial beta coeffs. during feedback in the conjunction ROIs for the cross-timepoint RSA

subInfo = readtable('Adult_Run_List.csv');
scriptPath = pwd; 

%% 2) Execute the roi searchlights

% IMPORTANT: specify rois here
roiPath = [pwd,'/ROI_combinedGLM','/conjunction_identity']; % change this for each model
roiList = dir(fullfile(roiPath,'*.nii'));

% IMPORTANT: specify your output directory here
output_folder = [pwd,'/Trialwise_Voxel_Betas/Feedback_combined_GLM'];

% IMPORTANT: specify which model you want to look at: 
% 1 = identity, 2 = qval
RDM_roi = 1; 

% define user options here 
userOptions = defineUserOptions_socialRisk(scriptPath);

for subIdx = 1:30

    % subjects that will be included in final analysis 
    if subInfo.exclusionCode(subIdx) == 0 
        
    % extract out the relevant onset info for both runs
    subNum = subInfo.subNum(subIdx); 
    
    % run for TG & SM
    for condition = 1:2

        % construct predictor RDMs 
        [pred_RDM,conditionVec] = construct_pred_RDMs_RSA_trialwise_feedback(subNum,condition,scriptPath,RDM_roi); % change this so you're only running relevant RDM
        extract_roi_coeffs_trialwise(subNum,scriptPath,pred_RDM,userOptions,condition,conditionVec,roiPath,roiList,output_folder,RDM_roi)
    end 
    
    elseif subInfo.exclusionCode(subIdx) == 1
        fprintf('Invalid Subject'); 
    end 
    
end

