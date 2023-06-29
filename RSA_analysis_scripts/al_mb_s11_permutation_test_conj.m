%% Run the non-parametric perumutation test on phase-combined results (feedback + choice)
% using SnPM package 

%% Step 1) import images and specify perumutations

subInfo = readtable('Adult_Run_List.csv');
scriptPath = pwd; 
combinedPath = [pwd,'/CombinedGLM_Feedback_3mm_Multimodel_Searchlight_9mm_1ss'];
secondLevelFolder = [pwd,'/Conj_CombinedGLM_Choice_3mm_Multimodel_Searchlight_9mm_1ss_Level2'];

method_type = {'b_est','corr_dist','corr_dist_ztrans','t_stat'};

% define searchlight radius
searchlight_radius = 9; 
step_size = 1;

% define the prediction RDM names here
condition = 1; 
subNum = 4002; % just grab from a random participant to get field names 
[pred_RDM,conditionVec] = construct_pred_RDMs_feedback(subNum,condition,scriptPath);  % make sure you're grabbing the right matrices 
modelNames = fieldnames(pred_RDM);

for RDM = 1:length(modelNames)

        for m = 1:length(method_type)

        method_str = char(method_type(m));
        
        % navigate to the relevant analysis directory 
        cur_folder = sprintf('%s/%s/%s/conj_test',secondLevelFolder, method_str,modelNames{RDM});
        
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
            subject_combined_folder = sprintf('/%d_searchlight_combined',subNum);
            outputPath = [combinedPath,subject_combined_folder];

            % find the relevant contrast
            conFilePath = sprintf('%s/s_sub-%d_%s_%s_combined.nii',...
                          outputPath,subNum,modelNames{RDM},method_str); 

            % add nifti image in
            scans{i} = conFilePath;

            else
            fprintf('Invalid Subject'); 
            end     

        end

        % will import the correct scans into group-level image 
        % important - use mask from choice phase for perm testing to get
        % conj results 
        matlabbatch{1}.spm.tools.snpm.des.OneSampT.P = scans';
        matlabbatch{1}.spm.tools.snpm.des.OneSampT.cov = struct('c', {}, 'cname', {});
        matlabbatch{1}.spm.tools.snpm.des.OneSampT.nPerm = 5000;
        matlabbatch{1}.spm.tools.snpm.des.OneSampT.vFWHM = [0 0 0];
        matlabbatch{1}.spm.tools.snpm.des.OneSampT.bVolm = 1;
        matlabbatch{1}.spm.tools.snpm.des.OneSampT.ST.ST_U = 0.0001;
        matlabbatch{1}.spm.tools.snpm.des.OneSampT.masking.tm.tm_none = 1;
        matlabbatch{1}.spm.tools.snpm.des.OneSampT.masking.im = 1;
        matlabbatch{1}.spm.tools.snpm.des.OneSampT.masking.em = {'/gpfs/data/ofeldman/alamba1/adult_social_risk/Final_pipeline_scripts/Final_RSA_analysis_scripts/Perm_Masks/choice_identity_tstat_perm_05_binary.nii'};
        matlabbatch{1}.spm.tools.snpm.des.OneSampT.globalc.g_omit = 1;
        matlabbatch{1}.spm.tools.snpm.des.OneSampT.globalm.gmsca.gmsca_no = 1;
        matlabbatch{1}.spm.tools.snpm.des.OneSampT.globalm.glonorm = 1;
        spm_jobman('run', matlabbatch); 

        end 
      
end 

cd(scriptPath);
clear all


%% Step 2) Compute cluster mass  

% you'll first need to create a list of .nii files in separate excel file 
scriptPath = pwd; 
combinedPath = [pwd,'/CombinedGLM_Feedback_3mm_Multimodel_Searchlight_9mm_1ss'];
secondLevelFolder = [pwd,'/Conj_CombinedGLM_Choice_3mm_Multimodel_Searchlight_9mm_1ss_Level2'];

method_type = {'b_est','corr_dist','corr_dist_ztrans','t_stat'};

% define searchlight radius
searchlight_radius = 9; 
step_size = 1;

% define the prediction RDM names here
condition = 1; 
subNum = 4002; % just grab from a random participant to get field names 
[pred_RDM,conditionVec] = construct_pred_RDMs_feedback(subNum,condition,scriptPath);  % make sure you're grabbing the right matrices 
modelNames = fieldnames(pred_RDM);

for RDM = 1:length(modelNames)

    for m = 1:length(method_type)

    method_str = char(method_type(m));
    
    % navigate to the relevant analysis directory 
    matlabbatch = {};
    cur_folder = sprintf('%s/%s/%s/conj_test',secondLevelFolder, method_str,modelNames{RDM});
    snpm_file = [cur_folder,'/SnPMcfg.mat'];
    matlabbatch{1}.spm.tools.snpm.cp.snpmcfg = {snpm_file};
    spm_jobman('run', matlabbatch); 

    end    
end 

%% Step 3) Inference: Note: this step needs to be done manually through the GUI 






