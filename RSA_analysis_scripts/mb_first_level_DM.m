%% Code constructs the first-level GLM design matrix for each subject 
% GLM includes choice and feedback for each trial, and with each partner
% type in the same GLM

function mb_first_level_DM(subNum,counterBalance,scriptPath)

%specify directory to store SPM.mat file (design matrix)
timePath = [scriptPath,'/Full_DM_Timings_Combined/'];
cd(timePath)

A1 = subNum;
str_subj = "%d_timings";
sub = sprintf(str_subj, A1); 
sub = char(sub);

str_tgDat = "TG_%d"; 
tgDat = sprintf(str_tgDat,A1);
tgDat = char(tgDat);

directoryPath = fullfile(timePath, sub, tgDat);
char(directoryPath);

spmPath = fullfile(timePath, sub);
char(spmPath);

matlabbatch = {};
matlabbatch{1}.spm.stats.fmri_spec.dir = {spmPath};
matlabbatch{1}.spm.stats.fmri_spec.timing.units = 'secs';
matlabbatch{1}.spm.stats.fmri_spec.timing.RT = 2;
matlabbatch{1}.spm.stats.fmri_spec.timing.fmri_t = 16;
matlabbatch{1}.spm.stats.fmri_spec.timing.fmri_t0 = 8;
%% TG run1 and run2 parameters are specified here 
%specify path to subject-specific NIFTI file 

niftyPath = [scriptPath,'/NIFTI_2mm/',int2str(subNum)];
runPath = dir(fullfile(niftyPath, '*ep2d_functional*'));
runPath  = struct2cell(runPath);

%% import TG session scans 
%need to factor in location of folder depending on counterbalanced order

if counterBalance == 1
    TG_run1Path = runPath(1,1);
    TG_run1Path  = cell2mat(TG_run1Path); 

    TG_run2Path = runPath(1,2);
    TG_run2Path  = cell2mat(TG_run2Path);

elseif counterBalance == 2
    TG_run1Path = runPath(1,3);
    TG_run1Path  = cell2mat(TG_run1Path); 

    TG_run2Path = runPath(1,4);
    TG_run2Path  = cell2mat(TG_run2Path); 
end

% import run 1 scans 
TG_run1_scans = {};
TG_run1_scannames = dir(fullfile(niftyPath,TG_run1Path, 'swra*.img'));

TG_R1_fullScanPath = fullfile(niftyPath,TG_run1Path);

    TG_run1_scans = [TG_run1_scans; ...
        reshape(strcat(TG_R1_fullScanPath, '/', {TG_run1_scannames.name}, ',1'), [], 1)];

% import run 2 scans 
TG_run2_scans = {};
TG_run2_scannames = dir(fullfile(niftyPath,TG_run2Path, 'swra*.img'));

TG_R2_fullScanPath = fullfile(niftyPath,TG_run2Path);

    TG_run2_scans = [TG_run2_scans; ...
        reshape(strcat(TG_R2_fullScanPath, '/', {TG_run2_scannames.name}, ',1'), [], 1)];

% combine scans 
TG_scans = [TG_run1_scans; TG_run2_scans];
    
matlabbatch{1}.spm.stats.fmri_spec.sess(1).scans = TG_scans;

%% specify TG condition parameters 

% load the subject timing parameters 
mcfile_str = "Cmbd_TG_Onsets_%d.mat"; 
mcPath = fullfile(directoryPath, char(sprintf(mcfile_str,A1)));
load(mcPath)

% restructure onset files
names = reshape(names,1,120);
durations = reshape(durations,1,120);
onsets = reshape(onsets,1,120);

for cond = 1:length(names)

matlabbatch{1}.spm.stats.fmri_spec.sess(1).cond(cond).name = names{cond};
matlabbatch{1}.spm.stats.fmri_spec.sess(1).cond(cond).onset = onsets{cond};
matlabbatch{1}.spm.stats.fmri_spec.sess(1).cond(cond).duration = durations{cond};
matlabbatch{1}.spm.stats.fmri_spec.sess(1).cond(cond).tmod = 0;

matlabbatch{1}.spm.stats.fmri_spec.sess(1).cond(cond).pmod = struct('name', {}, 'param', {}, 'poly', {});
matlabbatch{1}.spm.stats.fmri_spec.sess(1).cond(cond).orth = 1;

end 

%regressors 
matlabbatch{1}.spm.stats.fmri_spec.sess(1).regress = struct('name', {}, 'val', {});
regressor_str = "TG_combined_motionParams_%d.mat";
regressorPath = fullfile(directoryPath, char(sprintf(regressor_str,A1)));
matlabbatch{1}.spm.stats.fmri_spec.sess(1).multi_reg = {regressorPath};
matlabbatch{1}.spm.stats.fmri_spec.sess(1).hpf = 128;


%% import SM session scans 
% respecify timings path 
timePath = [scriptPath,'/Full_DM_Timings_Combined/'];
cd(timePath)

A1 = subNum;
str_subj = "%d_timings";
sub = sprintf(str_subj, A1); 
sub = char(sub);

str_smDat = "SM_%d"; 
smDat = sprintf(str_smDat,A1);
smDat = char(smDat);

directoryPath = fullfile(timePath, sub, smDat);
char(directoryPath);

niftyPath = [scriptPath,'/NIFTI_2mm/',int2str(subNum)];
runPath = dir(fullfile(niftyPath, '*ep2d_functional*'));
runPath  = struct2cell(runPath);

%need to factor in location of folder depending on counterbalanced order
if counterBalance == 1
    SM_run1Path = runPath(1,3);
    SM_run1Path  = cell2mat(SM_run1Path); 

    SM_run2Path = runPath(1,4);
    SM_run2Path  = cell2mat(SM_run2Path);
elseif counterBalance == 2
    SM_run1Path = runPath(1,1);
    SM_run1Path  = cell2mat(SM_run1Path); 

    SM_run2Path = runPath(1,2);
    SM_run2Path  = cell2mat(SM_run2Path); 
end

% import run 1 scans 
SM_run1_scans = {};
SM_run1_scannames = dir(fullfile(niftyPath,SM_run1Path, 'swra*.img'));

SM_R1_fullScanPath = fullfile(niftyPath,SM_run1Path);

    SM_run1_scans = [SM_run1_scans; ...
        reshape(strcat(SM_R1_fullScanPath, '/', {SM_run1_scannames.name}, ',1'), [], 1)];

% import run 2 scans 
SM_run2_scans = {};
SM_run2_scannames = dir(fullfile(niftyPath,SM_run2Path, 'swra*.img'));

SM_R2_fullScanPath = fullfile(niftyPath,SM_run2Path);

    SM_run2_scans = [SM_run2_scans; ...
        reshape(strcat(SM_R2_fullScanPath, '/', {SM_run2_scannames.name}, ',1'), [], 1)];

% combine scans 
SM_scans = [SM_run1_scans; SM_run2_scans];
    
matlabbatch{1}.spm.stats.fmri_spec.sess(2).scans = SM_scans;

%% specify SM condition parameters 

% load the subject timing parameters 
mcfile_str = "Cmbd_SM_Onsets_%d.mat"; 
mcPath = fullfile(directoryPath, char(sprintf(mcfile_str,A1)));
load(mcPath)

% restructure onset files
names = reshape(names,1,120);
durations = reshape(durations,1,120);
onsets = reshape(onsets,1,120);

for cond = 1:length(names) 

matlabbatch{1}.spm.stats.fmri_spec.sess(2).cond(cond).name = names{cond};
matlabbatch{1}.spm.stats.fmri_spec.sess(2).cond(cond).onset = onsets{cond};
matlabbatch{1}.spm.stats.fmri_spec.sess(2).cond(cond).duration = durations{cond};
matlabbatch{1}.spm.stats.fmri_spec.sess(2).cond(cond).tmod = 0;

matlabbatch{1}.spm.stats.fmri_spec.sess(2).cond(cond).pmod = struct('name', {}, 'param', {}, 'poly', {});
matlabbatch{1}.spm.stats.fmri_spec.sess(2).cond(cond).orth = 1;

end 

% Regressors 
matlabbatch{1}.spm.stats.fmri_spec.sess(2).regress = struct('name', {}, 'val', {});
regressor_str = "SM_combined_motionParams_%d.mat";
regressorPath = fullfile(directoryPath, char(sprintf(regressor_str,A1)));
matlabbatch{1}.spm.stats.fmri_spec.sess(2).multi_reg = {regressorPath};
matlabbatch{1}.spm.stats.fmri_spec.sess(2).hpf = 128;


%% Other 1st level GLM parameters. Mostly set to defaults, except time derivatives
matlabbatch{1}.spm.stats.fmri_spec.fact = struct('name', {}, 'levels', {});
matlabbatch{1}.spm.stats.fmri_spec.bases.hrf.derivs = [0 0];
matlabbatch{1}.spm.stats.fmri_spec.volt = 1;
matlabbatch{1}.spm.stats.fmri_spec.global = 'None';
matlabbatch{1}.spm.stats.fmri_spec.mthresh = 0.8;
matlabbatch{1}.spm.stats.fmri_spec.mask = {''};
matlabbatch{1}.spm.stats.fmri_spec.cvi = 'AR(1)';
    
%% Run the batch :)
spm_jobman('run', matlabbatch);    
cd(scriptPath)
