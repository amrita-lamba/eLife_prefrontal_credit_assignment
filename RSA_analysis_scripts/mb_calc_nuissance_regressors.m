%% pull out movement regressors and for each subject 
function mb_calc_nuissance_regressors(subNum,counterBalance,runxrun_correction,scriptPath)

niftyPath = [scriptPath,'/NIFTI_2mm/',int2str(subNum)];

if runxrun_correction == 0

    realignPath = dir(fullfile(niftyPath, '*ep2d_functional*'));
    realignPath = struct2cell(realignPath);
    realignPath = realignPath(1,1);
    realignPath = cell2mat(realignPath);

    realignFile = dir(fullfile(niftyPath, realignPath, '*rp*.txt'));
    realignFile = struct2cell(realignFile);
    realignFile = realignFile(1,1);
    realignFile = cell2mat(realignFile);

    % need to navigate to the folder to load the realign file 
    cd(fullfile(niftyPath,realignPath))
    realignParams = load(realignFile);

else 

    rundirs = dir(fullfile(niftyPath, '*ep2d_functional*'));
    rundirs = struct2cell(rundirs);
    rundirs = rundirs(1,:);

    cd(niftyPath)

    realignParams = [];

    for ri = 1:4
        cd(rundirs{ri})
        
        realignPath = dir(fullfile(niftyPath, '*ep2d_functional*'));
        realignPath = struct2cell(realignPath);
        realignPath = realignPath(1,1);
        realignPath = cell2mat(realignPath);

        %find the motion param file 
        realignFile = dir(fullfile(niftyPath, rundirs{ri}, '*rp*.txt'));
        realignFile = struct2cell(realignFile);
        realignFile = realignFile(1,1);
        realignFile = cell2mat(realignFile);
        runParams = load(realignFile);

        %concatonate motion params across runs into single arrays 
        realignParams = [realignParams; runParams]; 
        cd(niftyPath)
    end 

    % need to navigate to the folder to load the realign file 
    cd(fullfile(niftyPath,realignPath))
    
end 


%% calculate motion derivatives for Friston 24-parameter model
% can use for excessive motion but not used in manuscript 
% first 6 are raw motion params
% next 6 are 1 time point before
% last 12 are first 12 squared 

% populate with previous trials
% for deriv_col = 7:12
%     
%     raw_col = deriv_col -  6; 
%     
%     for img = 1:length(realignParams)
%         
%         if img == 1
%         realignParams(1,deriv_col) = 0; 
%         else
%         
%         prev_img = img - 1;     
%         realignParams(img,deriv_col) = realignParams(img,raw_col) - realignParams(prev_img,raw_col);     
%             
%         end 
% 
%     end 
% 
% end 
% 
% % calculate squares 
% for deriv_col = 13:24
%     
%     raw_col = deriv_col - 12; 
%     realignParams(:,deriv_col) = realignParams(:,raw_col).^2;
%    
% end 

%% add in framewise displacement binary
% load in FD file
fdFile = dir(fullfile(niftyPath, realignPath, 'FD_threshold.txt'));
fdFile = struct2cell(fdFile);
fdFile = fdFile(1,1);
fdFile = cell2mat(fdFile);
fdThreshold = load(fdFile);

% add in framewise displacement 
fdParams = dir(fullfile(niftyPath, realignPath, 'FD_params.txt'));
fdParams = struct2cell(fdParams);
fdParams = fdParams(1,1);
fdParams = cell2mat(fdParams);
FD = load(fdParams);
realignParams(:,7) = FD; 

%% add in CSF time course as regressor 

% niftyPath = fullfile('/Users/astrocyte/Dropbox (Brown)/AMRITA/Adolescent-risk Study/Adult_Subject_Data_Analyses/fmri_wholeBrain_ROI_pipeline_3/NIFTI/', int2str(subNum));
% cd(niftyPath); 
% 
% % read in CSF time course data
% csfParams = dir(fullfile(niftyPath, 'CSF_timecourse.txt'));
% csfParams = struct2cell(csfParams);
% csfParams = csfParams(1,1);
% csfParams = cell2mat(csfParams);
% csf = load(csfParams);
% realignParams(:,7) = csf; 

%% construct a nuissance regressor column for each high motion frame

% find high motion frames
highMotion = find(fdThreshold == 1);
rp_size = size(realignParams);

for highMoScan = 1:length(highMotion)

    regCol = highMoScan + rp_size(2); 
    nuissanceCol = zeros(length(fdThreshold),1);
    nuissanceCol(highMotion(highMoScan)) = 1; 
    realignParams(:,regCol) = nuissanceCol; 

end

%% construct additional nuissance regressors for first and last frames

% % find first and last frames
% endFrames = find(epThreshold == 1);
% rp_size = size(realignParams);
% 
% for ep = 1:length(endFrames)
% 
%     regCol = ep + rp_size(2); 
%     nuissanceCol = zeros(length(epThreshold),1);
%     nuissanceCol(endFrames(ep)) = 1; 
%     realignParams(:,regCol) = nuissanceCol; 
% 
% end


%% Set up timings directories

A1 = subNum;
str_subj = "%d_timings";
sub = sprintf(str_subj, A1);
sub = char(sub); 

str_tgDat = "TG_%d"; 
tgDat = sprintf(str_tgDat,A1);
tgDat = char(tgDat);

str_smDat = "SM_%d"; 
smDat = sprintf(str_smDat,A1);
smDat = char(smDat);

timePath = [scriptPath,'/Full_DM_Timings_Choice/'];

%% Save TG Run 1 realignment parameters  

% extract relevant motion params from rp. preprocessing doc 
if counterBalance == 1
TG_run1_motionParams = realignParams(1:246,:);
R1 =  TG_run1_motionParams;
% add run regressor 
runNum = repelem(1,246); 
R1 = [R1 runNum'];

elseif counterBalance == 2 
TG_run1_motionParams = realignParams(493:738,:);
R1 =  TG_run1_motionParams;  
% add run regressor 
runNum = repelem(3,246); 
R1 = [R1 runNum'];
end 

% navigate back to the subject folder for saving 
cd(fullfile(timePath, sub, tgDat));

% save movement params 
movement_file = "TG_run1_motionParams_%d";
save(sprintf(movement_file, A1),'R1');  

%% Save TG Run 2 realignment parameters

% extract relevant motion params from rp. preprocessing doc 
if counterBalance == 1
TG_run2_motionParams = realignParams(247:492,:);
R2 =  TG_run2_motionParams;
% add run regressor 
runNum = repelem(2,246); 
R2 = [R2 runNum'];

elseif counterBalance == 2 
TG_run1_motionParams = realignParams(739:984,:);
R2 =  TG_run1_motionParams;    
% add run regressor 
runNum = repelem(4,246); 
R2 = [R2 runNum'];
end 

% navigate back to the subject folder for saving 
cd(fullfile(timePath, sub, tgDat));

% save movement params 
movement_file = "TG_run2_motionParams_%d";
save(sprintf(movement_file, A1),'R2');  


%% Save SM Run 1 realignment params 

% extract relevant motion params 
if counterBalance == 1
SM_run1_motionParams = realignParams(493:738,:);
R1 =  SM_run1_motionParams;
% add run regressor 
runNum = repelem(3,246); 
R1 = [R1 runNum'];

elseif counterBalance == 2 
SM_run1_motionParams = realignParams(1:246,:);
R1 =  SM_run1_motionParams;
% add run regressor 
runNum = repelem(1,246); 
R1 = [R1 runNum'];
end 

% navigate back to the subject folder for saving 
cd(fullfile(timePath, sub, smDat));

% save movement params 
movement_file = "SM_run1_motionParams_%d";
save(sprintf(movement_file, A1),'R1');  

%% Save SM Run 2 realignment params

% extract relevant motion params 
if counterBalance == 1
SM_run2_motionParams = realignParams(739:984,:);
R2 =  SM_run2_motionParams;
% add run regressor 
runNum = repelem(4,246); 
R2 = [R2 runNum'];

elseif counterBalance == 2
SM_run2_motionParams = realignParams(247:492,:);
R2 =  SM_run2_motionParams;
runNum = repelem(2,246); 
R2 = [R2 runNum'];
end

% navigate back to the subject folder for saving 
cd(fullfile(timePath, sub, smDat));

% save movement params 
movement_file = "SM_run2_motionParams_%d";
save(sprintf(movement_file, A1),'R2');  
cd(scriptPath)


