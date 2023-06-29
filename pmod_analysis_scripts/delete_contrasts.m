%% this script will delete old contrats constrats and beta estimates 
%(will keep orignal NIFTY files)

for subNum = 4001:4030

% cd to the Timings folder
homePath = pwd; 
timePath = [pwd,'/Simple_DM_Timings'];

A1 = subNum;
str_subj = "%d_timings";
sub = sprintf(str_subj, A1);
sub = char(sub); 

%% this will loop through each subject's functional run and delete the 
% intermediate files
cd(fullfile(timePath, sub));
delete('beta*.nii');
delete('con*.nii');
delete('mask*.nii');
delete('Res*.nii');
delete('RPV*.nii');
delete('SPM*.mat');
delete('spmT*.nii');
  
cd(homePath)

end 


