%% Code for batch preprocessing data for univariate pmod analyses

% Takes the al_proproc_univ script (steps for realignment, slice time correction, 
% normalization, and smoothing) and runs for each subject 

scriptPath = pwd;
subInfo = readtable('Adult_Run_List.csv');

%% this section will preprocess the data for each subject 

for sub = 4001:4030
    al_proproc_univ(sub)
    cd(scriptPath)
end 

%% delete the intermediate preprocessing files 

cd(scriptPath)
for sub = 4001:4030
    delete_interm_files(sub)
    cd(scriptPath)
end 
