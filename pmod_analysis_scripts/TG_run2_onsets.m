%% This function pulls out the stimulus onset and run2_durations for each subject for TG
% to be passed through 1st level GLM specifcations 

function TG_run2_onsets(subNum,counterBalance,scriptPath)
% eventually this will be the input to a function 
% extract out the subject's TG data for each run seperately 

% this script will also organize all the feedback onsets and will calculate
% the prediction error values for all trials 
mlFile = readtable('adult_behavData_2ca2lr.csv');

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

cd(fullfile(timePath, sub, tgDat));
load('run2_tgData.mat')

% % % change to relevant runNum & condition 
sub_idx = find(mlFile.subjNum == subNum & mlFile.runNum == 2 & mlFile.condition == 1);
sub_data = mlFile(sub_idx,:); 

run2_names = {'all_decisions','missedTrial_decision',...
              'all_feedback','missedTrial_feedback'};
run2_onsets = {};
run2_durations = {}; 
run2_param_names = {'all_delta','missedTrial_delta'};
run2_param_vals = {};

faceOnsets = run2_tgData.faceOnset + 492;
faceDurations = run2_tgData.trialRT;
allChoices = run2_tgData.endowment;

feedbackOnsets = run2_tgData.feedbackOnset + 492; 
feedbackDurations = repelem(2,30);

% need to extract trials with missing values (i.e. subject did not respond quickly enough)
timP = find(allChoices == -1);
run2_tgData.run2_shuffledList.playerNumber(timP) = 1;

% collapse across player types 
timCM = find(allChoices > -1);
run2_tgData.run2_shuffledList.playerNumber(timCM) = 2;

% change value of missed trials to mean so they're effectively excluded
% get rid of all -1s and change to 2 for the run2_durations 
faceDurations(faceDurations == -1) = 2.00; 
delta = sub_data.delta'; 

% sort choices by player type 
for choiceType = 1:2
    
    idx = find(run2_tgData.run2_shuffledList.playerNumber == choiceType); 
   
    if choiceType == 2
       % choice onsets 
       pOnset = faceOnsets(idx);
       run2_onsets{1} = pOnset; 

       pDuration = faceDurations(idx);
       run2_durations{1} = pDuration;
       
       % feedback onsets 
       pFeedbackOnset = feedbackOnsets(idx);
       run2_onsets{3} = pFeedbackOnset; 

       pFeedbackDuration = feedbackDurations(idx);
       run2_durations{3} = pFeedbackDuration;

       pdelta = delta(idx);
       run2_param_vals{3} = pdelta;    
   
    else 
        
       pMOnset = faceOnsets(idx);
       run2_onsets{2} = pMOnset; 

       pMDuration = faceDurations(idx);
       run2_durations{2} = pMDuration;
       
       % feedback onsets 
       pMFeedbackOnset = feedbackOnsets(idx);
       run2_onsets{4} = pMFeedbackOnset; 

       pMFeedbackDuration = feedbackDurations(idx);
       run2_durations{4} = pMFeedbackDuration;

       pMdelta = delta(idx);
       run2_param_vals{4} = pMdelta; 
     
   
     end 
   
end
 
file_name = "MC_TG_run2_%d";
save(sprintf(file_name, A1),'run2_names','run2_durations','run2_onsets','run2_param_names','run2_param_vals');   
  
  
  