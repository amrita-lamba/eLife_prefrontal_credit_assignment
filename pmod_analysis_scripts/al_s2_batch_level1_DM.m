%% Generate a design matrix from raw data for each participant
% loops through data for each subject and extracts the onsets, durations and 
% movement regressors for each subject and saves them in the
% appropriate format within the specific subject's timing folder 

subInfo = readtable('Adult_Run_List.csv');
scriptPath = pwd; 

%%  extract the relevant onsets and motion params 
% construct the design matrix per subject

for subIdx = 1:30
    
    % subjects that will be included in final analysis 
    if subInfo.exclusionCode(subIdx) == 0 
        
    % extract out the relevant onset info for both runs
    subNum = subInfo.subNum(subIdx); 
    counterBalance = subInfo.counterBalance(subIdx);
    runxrun_correction = subInfo.runxrun_correction(subIdx);
    
    % create motion regressors 
    calc_FD(subNum,runxrun_correction, scriptPath)
    calc_nuissance_regressors(subNum,counterBalance,runxrun_correction,scriptPath);
    cd(scriptPath);
    
    % create onset files for each TG run here 
    TG_run1_onsets(subNum,counterBalance,scriptPath);
    cd(scriptPath);
    TG_run2_onsets(subNum,counterBalance,scriptPath);
    cd(scriptPath);
    
    % create onset files for each SM run here 
    SM_run1_onsets(subNum,counterBalance,scriptPath);
    cd(scriptPath);
    SM_run2_onsets(subNum,counterBalance,scriptPath);
    cd(scriptPath);
    
    % combine onset files across runs here 
    combine_run_onsets(subNum,scriptPath); 
    cd(scriptPath);
    
    elseif subInfo.exclusionCode(subIdx) == 1
        fprintf('Invalid Subject'); 
    end 
end 

%% Constructs the first-level DM 

for subIdx = 1:30
    
    if subInfo.exclusionCode(subIdx) == 0 
        
    % extract out the relevant onset info for both runs
    subNum = subInfo.subNum(subIdx); 
    counterBalance = subInfo.counterBalance(subIdx);
    first_level_DM(subNum,counterBalance,scriptPath);  
    
    elseif subInfo.exclusionCode(subIdx) == 1
        fprintf('Invalid Subject'); 
    end 
    cd(scriptPath);  

end 

%% Estimate the design matrix 

for subIdx =  1:30
    %subjects that will be included in final analysis 
    if subInfo.exclusionCode(subIdx) == 0 
       
    subNum = subInfo.subNum(subIdx); 
    estimate_DM(subNum,scriptPath);
    
    elseif subInfo.exclusionCode(subIdx) == 1
        fprintf('Invalid Subject'); 
    end 
    cd(scriptPath);  
end 

