%% This script will estimate the first level GLM for each subject 
function mb_estimate_DM(subNum,scriptPath)
%% Path to relevant folder 

timePath = [scriptPath,'/Full_DM_Timings_Combined/'];

A1 = subNum;
str_subj = "%d_timings";
sub = sprintf(str_subj, A1);
sub = char(sub);
dmPath = (fullfile(timePath, sub, 'SPM.mat'));

%% specify estimation parameters 
matlabbatch = {};
matlabbatch{1}.spm.stats.fmri_est.spmmat = {dmPath};
matlabbatch{1}.spm.stats.fmri_est.write_residuals = 0;
matlabbatch{1}.spm.stats.fmri_est.method.Classical = 1;

%% Run the batch :)
spm_jobman('run', matlabbatch); 
cd(scriptPath)

