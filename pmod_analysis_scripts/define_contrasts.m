%% Defines the contrasts (basic t test) for first level analyses
% set 1 & 2 are basic contrasts during choice & feedback phase 

function define_contrasts(subNum,scriptPath)

% import the contrast codes 
con_code_ML = readtable('design_matrix_code.csv');

% update the relevant columns here if .csv file changes 
TG_code_columns = 4:8; 
SM_code_columns = 9:13; 

% cd to the Timings folder
timePath = fullfile(scriptPath,'/Simple_DM_Timings/');
cd(timePath)

A1 = subNum;
str_subj = "%d_timings";
sub = sprintf(str_subj, A1);
sub = char(sub); 

str_tgDat = "TG_%d"; 
tgDat = sprintf(str_tgDat,A1);
tgDat = char(tgDat);

% cd to .SPM folder
cd(fullfile(timePath, sub));

matlabbatch = [];
matlabbatch{1}.spm.stats.con.spmmat =  cellstr(strcat(fullfile(timePath, sub),'/SPM.mat'));

% pull out the subject's motion parameters file
cd(fullfile(timePath, sub, tgDat));
motionFile = dir(fullfile(timePath, sub, tgDat, 'TG_combined*.mat'));
motionFile = struct2cell(motionFile);
motionFile = motionFile(1,1);
motionFile = cell2mat(motionFile);
load(motionFile);

% specific regressor vector
R_size = size(R); 
regressor_vec = zeros(1,R_size(2));

%% Construct the contrasts from the design matrix code list 

for conNum = 1:10
             
    % specify TG contrast vector 
    TG_vec = con_code_ML(conNum,TG_code_columns); 
    TG_vec = table2array(TG_vec);
    
    % specify SM contrast vector 
    SM_vec = con_code_ML(conNum,SM_code_columns); 
    SM_vec = table2array(SM_vec);
    conVec = [TG_vec, regressor_vec, SM_vec, regressor_vec]; 
   
    % specify contrast name 
    conName = con_code_ML.conName(conNum);
    conName = char(conName);

    % create the batch 
    matlabbatch{1}.spm.stats.con.consess{conNum}.tcon.name = conName;
    matlabbatch{1}.spm.stats.con.consess{conNum}.tcon.convec = conVec;
    matlabbatch{1}.spm.stats.con.consess{conNum}.tcon.sessrep = 'none';

end 

%% Run the contrasts 
spm_jobman('run',matlabbatch); 

