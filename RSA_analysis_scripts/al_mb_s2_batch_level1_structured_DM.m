%% Generate a design matrix from raw data for each subject 
%loops through data for each subject and extracts the onsets, durations and 
%movement regressors for each subject and saves them in the
%appropriate format within the specific subject's timing folder 

%read in subject info table with counter-balance and exclusion information
subInfo = readtable('Adult_Run_List.csv');

%scriptPath with the relevant functions needed to construct the design matrix
%make sure you loop back here in order to index the correct script s
scriptPath = pwd; 

%% Get motion, onsets, and duration params for each subject
for subIdx = 1:30
    % subjects that will be included in final analysis 
    if subInfo.exclusionCode(subIdx) == 0 
        
    % get relevant subject info
    subNum = subInfo.subNum(subIdx); 
    counterBalance = subInfo.counterBalance(subIdx);
    runxrun_correction = subInfo.runxrun_correction(subIdx);
    
    % create motion regressors 
    mb_calc_FD(subNum,runxrun_correction,scriptPath);
    mb_calc_nuissance_regressors(subNum,counterBalance,runxrun_correction,scriptPath);
    
    % combine motion regressors
    mb_combine_motion_regressors(subNum,scriptPath);
    
    % create onset files for each TG run here 
    mb_TG_onsets_combined(subNum,counterBalance,scriptPath);
    mb_SM_onsets_combined(subNum,counterBalance,scriptPath);
    
    elseif subInfo.exclusionCode(subIdx) == 1
        fprintf('Invalid Subject'); 
    end 
end 

%% Construct the first-level DM 

for subIdx = 1:30
    
    if subInfo.exclusionCode(subIdx) == 0 
        
    % extract out the relevant onset info for both runs
    subNum = subInfo.subNum(subIdx); 
    counterBalance = subInfo.counterBalance(subIdx);
    mb_first_level_DM(subNum,counterBalance,scriptPath);  
    
    elseif subInfo.exclusionCode(subIdx) == 1
        fprintf('Invalid Subject'); 
    end  

end 

%% Estimate the design matrix 

for subIdx =  1:30
    %subjects that will be included in final analysis 
    if subInfo.exclusionCode(subIdx) == 0 
       
    subNum = subInfo.subNum(subIdx); 
    mb_estimate_DM(subNum,scriptPath);
    
    elseif subInfo.exclusionCode(subIdx) == 1
        fprintf('Invalid Subject'); 
    end 
end 

