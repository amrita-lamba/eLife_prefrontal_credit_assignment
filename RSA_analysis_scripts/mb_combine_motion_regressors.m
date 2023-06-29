function mb_combine_motion_regressors(subNum,scriptPath)
%% Get TG motion params

timePath = [scriptPath,'/Full_DM_Timings_Combined/'];
cd(timePath)

A1 = subNum;
str_subj = "%d_timings";
sub = sprintf(str_subj, A1);
sub = char(sub); 

str_tgDat = "TG_%d"; 
tgDat = sprintf(str_tgDat,A1);
tgDat = char(tgDat);
cd(fullfile(timePath, sub, tgDat));

% combine movement parameters into a single file 
TG_run1_motion = "TG_run1_motionParams_%d";
TG_run1_motion_file = sprintf(TG_run1_motion,A1);
load(TG_run1_motion_file);

TG_run2_motion = "TG_run2_motionParams_%d";
TG_run2_motion_file = sprintf(TG_run2_motion,A1);
load(TG_run2_motion_file);

R = [R1;R2];

% save combined movement params 
movement_file = "TG_combined_motionParams_%d";
save(sprintf(movement_file, A1),'R');  

%% repeat for SM

cd(timePath)

A1 = subNum;
str_subj = "%d_timings";
sub = sprintf(str_subj, A1);
sub = char(sub); 

str_smDat = "SM_%d"; 
smDat = sprintf(str_smDat,A1);
smDat = char(smDat);
cd(fullfile(timePath, sub, smDat));

% combine movement parameters into a single file 
SM_run1_motion = "SM_run1_motionParams_%d";
SM_run1_motion_file = sprintf(SM_run1_motion,A1);
load(SM_run1_motion_file);

SM_run2_motion = "SM_run2_motionParams_%d";
SM_run2_motion_file = sprintf(SM_run2_motion,A1);
load(SM_run2_motion_file);

R = [R1;R2];

% save combined movement params 
movement_file = "SM_combined_motionParams_%d";
save(sprintf(movement_file, A1),'R');  
cd(scriptPath)

end 