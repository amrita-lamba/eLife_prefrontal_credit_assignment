%% Constructs the relevant level-1 univariate contrasts for each subject 
 
subInfo = readtable('Adult_Run_List.csv');
scriptPath = pwd; 

%% this section will run the basic set of simplified contrasts

for subIdx = 1:30
    if subInfo.exclusionCode(subIdx) == 0 
        subNum = subInfo.subNum(subIdx); 
        define_contrasts(subNum,scriptPath);
        cd(scriptPath);
    else 
        fprintf('Invalid Subject'); 
    end   
end 

