%% Run the non-parametric perumutation test on combined TG + SM delta image


%% Step 1) Import images and specify perumutations

subInfo = readtable('Adult_Run_List.csv');
scriptPath = pwd; 
combinedPath = [pwd,'/Simple_DM_Timings'];
secondLevelFolder = [pwd,'/Simple_DM_Timings_Level_2'];

% navigate to the relevant analysis directory 
cur_folder = sprintf('%s/perm_test',secondLevelFolder);

if ~exist(cur_folder)
    mkdir(cur_folder)
end

% intialize empty scan vector
scans = {};
matlabbatch = {};
matlabbatch{1}.spm.tools.snpm.des.OneSampT.DesignName = 'MultiSub: One Sample T test on diffs/contrasts';
matlabbatch{1}.spm.tools.snpm.des.OneSampT.DesignFile = 'snpm_bch_ui_OneSampT';
matlabbatch{1}.spm.tools.snpm.des.OneSampT.dir = {cur_folder};

%find that subject's .con file
i = 0; 
for sub = 1:30 

    if subInfo.exclusionCode(sub) == 0 

    % increment counter 
    i = i + 1;   
    subNum = subInfo.subNum(sub);

    % navigate to subject's director & select contrast
    subject_combined_folder = sprintf('/%d_timings',subNum);
    outputPath = [combinedPath,subject_combined_folder];

    % find the relevant contrast
    conFilePath = sprintf('%s/con_0010.nii',outputPath); % just running on delta image 

    % add nifti image in
    scans{i} = conFilePath;

    else
    fprintf('Invalid Subject'); 
    end     

end

% will import the correct scans into group-level image 
matlabbatch{1}.spm.tools.snpm.des.OneSampT.P = scans';
matlabbatch{1}.spm.tools.snpm.des.OneSampT.cov = struct('c', {}, 'cname', {});
matlabbatch{1}.spm.tools.snpm.des.OneSampT.nPerm = 5000;
matlabbatch{1}.spm.tools.snpm.des.OneSampT.vFWHM = [0 0 0];
matlabbatch{1}.spm.tools.snpm.des.OneSampT.bVolm = 1;
matlabbatch{1}.spm.tools.snpm.des.OneSampT.ST.ST_U = 0.0001;
matlabbatch{1}.spm.tools.snpm.des.OneSampT.masking.tm.tm_none = 1;
matlabbatch{1}.spm.tools.snpm.des.OneSampT.masking.im = 1;
matlabbatch{1}.spm.tools.snpm.des.OneSampT.masking.em = {scriptPath,'Subject_Masks/group_whole_brain_binary.nii'};
matlabbatch{1}.spm.tools.snpm.des.OneSampT.globalc.g_omit = 1;
matlabbatch{1}.spm.tools.snpm.des.OneSampT.globalm.gmsca.gmsca_no = 1;
matlabbatch{1}.spm.tools.snpm.des.OneSampT.globalm.glonorm = 1;
spm_jobman('run', matlabbatch); 

cd(scriptPath);
clear all


%% Step 2) Compute cluster mass 

subInfo = readtable('Adult_Run_List.csv');
scriptPath = pwd; 
combinedPath = [pwd,'/Simple_DM_Timings'];
secondLevelFolder = [pwd,'/Simple_DM_Timings_Level_2'];

% navigate to the relevant analysis directory 
matlabbatch = {};
cur_folder = sprintf('%s/perm_test',secondLevelFolder);
snpm_file = [cur_folder,'/SnPMcfg.mat'];
matlabbatch{1}.spm.tools.snpm.cp.snpmcfg = {snpm_file};
spm_jobman('run', matlabbatch); 

%% Step 3) Inference: Note: this step needs to be done manually through the GUI 






