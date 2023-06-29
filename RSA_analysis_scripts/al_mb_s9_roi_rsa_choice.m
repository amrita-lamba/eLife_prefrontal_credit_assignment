%% extract coefficents for each participant in specified ROIs during choice phase 

%read subject info table with counter-balance and exclusion data 
subInfo = readtable('Adult_Run_List.csv');
scriptPath = pwd; 

%% 2) Execute the roi searchlights

% IMPORTANT: specify rois here
roiPath = [pwd,'/ROI_combinedGLM','/choice_identity']; % change this for each model
roiList = dir(fullfile(roiPath,'*.nii'));

% IMPORTANT: specify your output directory here
output_folder = [pwd,'/CombinedGLM_Choice_3mm_Multimodel_Searchlight_9mm_1ss'];

% IMPORTANT: specify which model you want to look at: 
% 1 = identity, 2 = qval
RDM_roi = 1; 

% define user options here 
userOptions = defineUserOptions_socialRisk(scriptPath);

for subIdx = 1:30

    % subjects that will be included in final analysis 
    if subInfo.exclusionCode(subIdx) == 0 

    subNum = subInfo.subNum(subIdx); 

    % run for TG & SM
    for condition = 1:2

        [pred_RDM,conditionVec] = construct_pred_RDMs_choice(subNum,condition,scriptPath); 
        extract_roi_coeffs(subNum,scriptPath,pred_RDM,userOptions,condition,conditionVec,roiPath,roiList,output_folder,RDM_roi)
        
    end 
    
    elseif subInfo.exclusionCode(subIdx) == 1
        fprintf('Invalid Subject'); 
    end 
    
end






