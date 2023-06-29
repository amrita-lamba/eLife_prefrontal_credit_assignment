% Code for batch preprocessing of subject data 
% Takes the al_proproc script (steps for slice time correction, 
% realignment,normalization, and smoothing) and runs for each subject folder

scriptPath = pwd; 
subInfo = readtable('Adult_Run_List.csv');

%% Preprocess the data for each subject 

%change for relevant subjects 
for sub = 4001:4030
    al_proproc_RSA(sub,scriptPath)
    cd(scriptPath)
end 

%% Delete all the intermediate preprocessing files 

cd(scriptPath)
for sub = 4001:4030
    delete_interm_files(sub,scriptPath)
    cd(scriptPath)
end 

%% Note: The following steps are only necessary if motion/distortion is an issue

%%  Generate probalistic maps of CSF from anatomical 

cd(scriptPath)
for sub = 4001:4030
    segment_CSF(sub,scriptPath)
    cd(scriptPath)
end 

%% Create ROIs and pull out time course for each subject's CSF

cd(scriptPath)
for sub = 4001:4030
    timecourse_CSF(sub,scriptPath)
    cd(scriptPath)
end 


