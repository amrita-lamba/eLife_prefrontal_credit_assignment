%% Take the TG_run1, TG_run2, SM_run1 & SM_run2 files 
% and will combine runs together for each condition. Will result in an
% onset file for each condition (TG & SM) to input in the Design Matrix 

function combine_run_onsets(subNum,scriptPath) 
%% Combine TG Onset files  

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

% read in TG run1 parameters 
cd(fullfile(timePath, sub, tgDat));
TG_run1 = "MC_TG_Run1_%d";
TG_run1_file = sprintf(TG_run1,A1);
load(TG_run1_file);

% read in TG run2 parameters 
TG_run2 = "MC_TG_Run2_%d";
TG_run2_file = sprintf(TG_run2,A1);
load(TG_run2_file);

%%%%%%%% MEAN CENTER DATA %%%%%%%%%%
% get the mean of all choices, excluding missed trials 
TG_run1_vals = cell2mat(run1_param_vals); 
TG_run2_vals = cell2mat(run2_param_vals); 

TG_run1_returned = TG_run1_vals(1:30);
TG_run2_returned = TG_run2_vals(1:30);
allReturned = [TG_run1_returned, TG_run2_returned];

% find mean of each parametric modulator 
TG_mean_returned = mean(allReturned(allReturned ~= -1));

% mean center all choices excluding the missed trials 
% missed trials should always result in a coded value of 0
feedback_stDev = std(allReturned(allReturned ~= TG_mean_returned));

% z score run 1 feedback 
data = run1_param_vals{3};
run1_param_vals{3} = (data - TG_mean_returned)/feedback_stDev;

% % z score run 2 feedback 
data = run2_param_vals{3};
run2_param_vals{3} = (data - TG_mean_returned)/feedback_stDev;

% save all parameters into an array
durations = {};
names = {};
onsets = {};
param_names = {};
param_vals = {};

% aggregate values across runs 
for i = 1:4

durations{i} = [run1_durations{i},run2_durations{i}]; 
onsets{i} = [run1_onsets{i},run2_onsets{i}];  
names{i} = run1_names{i};

end 

param_names{1} = [run1_param_names{1}]; 
param_vals{1} = [run1_param_vals{3},run2_param_vals{3}]; 

TG_file_name = "Cmbd_TG_Onsets_%d";
save(sprintf(TG_file_name, A1),'names','durations','onsets','param_names','param_vals');   

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

% cd to the Timings folder
timePath = fullfile(scriptPath,'/Simple_DM_Timings/');
cd(timePath)

A1 = subNum;
str_subj = "%d_timings";
sub = sprintf(str_subj, A1);
sub = char(sub); 
str_smDat = "SM_%d"; 
smDat = sprintf(str_smDat,A1);
smDat = char(smDat);

% read in TG run1 parameters 
cd(fullfile(timePath, sub, smDat));
SM_run1 = "MC_SM_Run1_%d";
SM_run1_file = sprintf(SM_run1,A1);
load(SM_run1_file);

% read in TG run2 parameters 
SM_run2 = "MC_SM_Run2_%d";
SM_run2_file = sprintf(SM_run2,A1);
load(SM_run2_file);

%%%%%%%% MEAN CENTER DATA %%%%%%%%%%
% get the mean of all choices, excluding missed trials 
SM_run1_vals = cell2mat(run1_param_vals); 
SM_run2_vals = cell2mat(run2_param_vals); 

SM_run1_returned = SM_run1_vals(1:30);
SM_run2_returned = SM_run2_vals(1:30);
allReturned = [SM_run1_returned, SM_run2_returned];

% find mean of each parametric modulator 
SM_mean_returned = mean(allReturned(allReturned ~= -1));
allReturned(allReturned == -1) = SM_mean_returned; 

% mean center all choices excluding the missed trials 
% missed trials should always result in a coded value of 0
feedback_stDev = std(allReturned(allReturned ~= SM_mean_returned));

% z score run 1 feedback 
data = run1_param_vals{3};
run1_param_vals{3} = (data - SM_mean_returned)/feedback_stDev;

% % z score run 2 feedback 
data = run2_param_vals{3};
run2_param_vals{3} = (data - SM_mean_returned)/feedback_stDev;

% save all parameters into an array
durations = {};
names = {};
onsets = {};
param_names = {};
param_vals = {};

% aggregate values across runs 
for i = 1:4

durations{i} = [run1_durations{i},run2_durations{i}]; 
onsets{i} = [run1_onsets{i},run2_onsets{i}];  
names{i} = run1_names{i};

end 

param_names{1} = [run1_param_names{1}]; 
param_vals{1} = [run1_param_vals{3},run2_param_vals{3}]; 

SM_file_name = "Cmbd_SM_Onsets_%d";
save(sprintf(SM_file_name, A1),'names','durations','onsets','param_names','param_vals');   

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

