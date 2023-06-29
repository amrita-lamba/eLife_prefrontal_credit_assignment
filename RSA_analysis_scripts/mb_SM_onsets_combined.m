%% This function pulls out the stimulus onset and run1_durations for each subject for sm
% to be passed through 1st level GLM specifcations 
function mb_SM_onsets_combined(subNum,counterBalance,scriptPath)

mlFile = readtable('adult_behavData_2ca2lr.csv');

choice_conditionNames = readtable('betaSeries_fullConditionNames_choice.xlsx'); 
feedback_conditionNames = readtable('betaSeries_fullConditionNames_feedback.xlsx'); 

% cd to the Timings folder
timePath = [scriptPath,'/Full_DM_Timings_Combined/'];
cd(timePath)

A1 = subNum;
str_subj = "%d_timings";
sub = sprintf(str_subj, A1);
sub = char(sub); 

str_smDat = "SM_%d"; 
smDat = sprintf(str_smDat,A1);
smDat = char(smDat);

cd(fullfile(timePath, sub, smDat));
load('run1_smData.mat')
load('run2_smData.mat')

% % % change to relevant runNum & condition 
sub_idx = find(mlFile.subjNum == subNum & mlFile.condition == 2);
sub_data = mlFile(sub_idx,:); 

% get choice condition names, change for each run 
low_choice_names = choice_conditionNames.low_choice;
high_choice_names = choice_conditionNames.high_choice;
neutral_choice_names = choice_conditionNames.neutral_choice;
random_choice_names = choice_conditionNames.random_choice;

low_feedback_names = feedback_conditionNames.low_feedback;
high_feedback_names = feedback_conditionNames.high_feedback;
neutral_feedback_names = feedback_conditionNames.neutral_feedback;
random_feedback_names = feedback_conditionNames.random_feedback;
 
% aggregate
choice_names = [low_choice_names, high_choice_names, neutral_choice_names, random_choice_names]';
choice_onsets = {};
choice_durations = {}; 

feedback_names = [low_feedback_names, high_feedback_names, neutral_feedback_names, random_feedback_names]';
feedback_onsets = {};
feedback_durations = {}; 

% pull feedback onsets and durations 
run1_faceOnsets = run1_smData.faceOnset;
run1_faceDurations = run1_smData.trialRT;
run2_faceOnsets = run2_smData.faceOnset + 492;
run2_faceDurations = run2_smData.trialRT;

% pull feedback onsets and durations 
run1_feedbackOnsets = run1_smData.feedbackOnset; 
run2_feedbackOnsets = run2_smData.feedbackOnset + 492; 
feedbackDurations = repelem(2,60);

% combine run 1 & run 2 onsets
faceOnsets = [run1_faceOnsets, run2_faceOnsets];
faceDurations = [run1_faceDurations, run2_faceDurations]; 
feedbackOnsets = [run1_feedbackOnsets, run2_feedbackOnsets];

for trial = 1:60

binned_roundNum = sub_data.binned_trial_num(trial); 
playerNum = sub_data.playerNum(trial); 

    if playerNum == 1

           % choice onsets  
           choice_onsets{1,binned_roundNum} = faceOnsets(trial); 
           choice_durations{1,binned_roundNum} = faceDurations(trial);

           % feedback onsets 
           feedback_onsets{1,binned_roundNum} = feedbackOnsets(trial); 
           feedback_durations{1,binned_roundNum} = feedbackDurations(trial);  

    elseif playerNum == 2

           % choice onsets  
           choice_onsets{2,binned_roundNum} = faceOnsets(trial); 
           choice_durations{2,binned_roundNum} = faceDurations(trial);

           % feedback onsets 
           feedback_onsets{2,binned_roundNum} = feedbackOnsets(trial); 
           feedback_durations{2,binned_roundNum} = feedbackDurations(trial);  

    elseif playerNum == 3

           % choice onsets  
           choice_onsets{3,binned_roundNum} = faceOnsets(trial); 
           choice_durations{3,binned_roundNum} = faceDurations(trial);

           % feedback onsets 
           feedback_onsets{3,binned_roundNum} = feedbackOnsets(trial); 
           feedback_durations{3,binned_roundNum} = feedbackDurations(trial);  

    elseif playerNum == 4

           % choice onsets  
           choice_onsets{4,binned_roundNum} = faceOnsets(trial); 
           choice_durations{4,binned_roundNum} = faceDurations(trial);

           % feedback onsets 
           feedback_onsets{4,binned_roundNum} = feedbackOnsets(trial); 
           feedback_durations{4,binned_roundNum} = feedbackDurations(trial);  

    end 
end

% combine in interleaving order 
names = {}; 
onsets = {};
durations = {};

for i = 1:15   
    
    col = {choice_names{:,i};feedback_names{:,i}}';
    names = [names, col];
    
    vals = {choice_onsets{:,i};feedback_onsets{:,i}}';
    onsets = [onsets, vals];
   
    durs = {choice_durations{:,i};feedback_durations{:,i}}';
    durations = [durations, durs];

end


file_name = "Cmbd_SM_Onsets_%d";
save(sprintf(file_name, A1),'names','durations','onsets');   
cd(scriptPath)  
  