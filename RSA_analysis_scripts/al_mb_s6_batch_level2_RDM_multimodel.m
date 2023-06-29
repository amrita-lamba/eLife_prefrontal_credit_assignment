%%  Create group-level brain maps for RSA analysis

subInfo = readtable('Adult_Run_List.csv');
scriptPath = pwd; 

% define searchlight radius
searchlight_radius = 9; 
step_size = 1;

% define the prediction RDM names here
condition = 1; 
subNum = 4002; % just grab from a random participant to get field names 
[pred_RDM,conditionVec] = construct_pred_RDMs_choice(subNum,condition,scriptPath);  % make sure you're grabbing the right matrices 
modelNames = fieldnames(pred_RDM);

method_type = {'b_est','corr_dist','corr_dist_ztrans','t_stat'};

%% Step 1) Create level 2 .SPM files 

for condition = 1:2

    for RDM = 1:length(modelNames)
        
        if condition == 1
            condition_str = 'TG';
        else
            condition_str = 'SM';
        end 
        
        % loop over each distance method 
        for m = 1:length(method_type)

            method_str = char(method_type(m));

            % navigate to the relevant analysis directory 
            secondLevelFolder = sprintf('CombinedGLM_Choice_3mm_Multimodel_Searchlight_%dmm_%dss_Level2/%s/%s_%s',...
                searchlight_radius,step_size,method_str,condition_str,modelNames{RDM});
            cur_folder = fullfile(scriptPath, secondLevelFolder);    
            matlabbatch{1}.spm.stats.factorial_design.dir = {cur_folder};

            % intialize empty scan vector
            scans = {};

            %find that subject's .con file
            i = 0; 

            for sub = 1:30  

                if subInfo.exclusionCode(sub) == 0 

                % increment counter 
                i = i + 1;   
                subNum = subInfo.subNum(sub);

                % navigate to subject's director & select contrast
                niftiPath = sprintf('%s/CombinedGLM_Choice_3mm_Multimodel_Searchlight_%dmm_%dss/%d_searchlight_condition_%d',...
                            scriptPath,searchlight_radius,step_size,subNum,condition);

                conFile = dir(fullfile(niftiPath, 's_sub*.nii')); % smoothed/unsmoothed 
                conFile = struct2cell(conFile);
                conFile = conFile(1,1:length(modelNames));

                % find the relevant contrast
                conFilePath = sprintf('%s/s_sub-%d_%s_%s.nii',...
                              niftiPath,subNum,modelNames{RDM},method_str); % smoothed/unsmoothed

                % add nifti image in
                scans{i} = conFilePath;

                else
                fprintf('Invalid Subject'); 
                end     

            end

            % will import the correct scans into group-level image 
            matlabbatch{1}.spm.stats.factorial_design.des.t1.scans = scans';
            matlabbatch{1}.spm.stats.factorial_design.cov = struct('c', {}, 'cname', {}, 'iCFI', {}, 'iCC', {});
            matlabbatch{1}.spm.stats.factorial_design.multi_cov = struct('files', {}, 'iCFI', {}, 'iCC', {});
            matlabbatch{1}.spm.stats.factorial_design.masking.tm.tm_none = 1;
            matlabbatch{1}.spm.stats.factorial_design.masking.im = 1;
            matlabbatch{1}.spm.stats.factorial_design.masking.em = {'/gpfs/data/ofeldman/alamba1/adult_social_risk/fmri_modelBased_RSA/Subject_Masks_OLD/group_whole_brain_binary.nii,1'};
            matlabbatch{1}.spm.stats.factorial_design.globalc.g_omit = 1;
            matlabbatch{1}.spm.stats.factorial_design.globalm.gmsca.gmsca_no = 1;
            matlabbatch{1}.spm.stats.factorial_design.globalm.glonorm = 1;
            spm_jobman('run', matlabbatch); 
        
        end 
     end 
end 

cd(scriptPath);
clear all


%% Step 2) Estimate the model

scriptPath = pwd; 

% define searchlight radius
searchlight_radius = 9; 
step_size = 1;

% define the prediction RDM names here
condition = 1; 
subNum = 4002; % just grab from a random participant to get field names 
[pred_RDM,conditionVec] = construct_pred_RDMs_choice(subNum,condition,scriptPath);
modelNames = fieldnames(pred_RDM);

method_type = {'b_est','corr_dist','corr_dist_ztrans','t_stat'};

for condition = 1:2
    
    for RDM = 1:length(modelNames)
        
        if condition == 1
            condition_str = 'TG';
        else
            condition_str = 'SM';
        end
        
        % loop over each distance method
        for m = 1:length(method_type)

            method_str = char(method_type(m));

            % navigate to the relevant analysis directory 
            secondLevelFolder = sprintf('/CombinedGLM_Choice_3mm_Multimodel_Searchlight_%dmm_%dss_Level2/%s/%s_%s',...
                searchlight_radius,step_size,method_str,condition_str,modelNames{RDM});
            cur_folder = fullfile(scriptPath, secondLevelFolder);    

            % navigate to the relevant analysis directory 
            groupImFile = fullfile(cur_folder, 'SPM.mat');
            matlabbatch{1}.spm.stats.fmri_est.spmmat = {groupImFile};
            matlabbatch{1}.spm.stats.fmri_est.write_residuals = 0;
            matlabbatch{1}.spm.stats.fmri_est.method.Classical = 1;
            spm_jobman('run', matlabbatch); 
        
        end 
    end 
end 

cd(scriptPath);
clear all 

%% Step 3) T-test against 0

scriptPath = pwd; 

% define searchlight radius
searchlight_radius = 9; 
step_size = 1;

% define the prediction RDM names here
condition = 1; 
subNum = 4002; % just grab from a random participant to get field names 
[pred_RDM,conditionVec] = construct_pred_RDMs_choice(subNum,condition,scriptPath);
modelNames = fieldnames(pred_RDM);

method_type = {'b_est','corr_dist','corr_dist_ztrans','t_stat'};

for condition = 1:2
    
    for RDM = 1:length(modelNames)
        
        if condition == 1
            condition_str = 'TG';
        else
            condition_str = 'SM';
        end 
        
        % loop over each distance method 
        % loop over each distance method
        for m = 1:length(method_type)

            method_str = char(method_type(m));

            % navigate to the relevant analysis directory 
            secondLevelFolder = sprintf('/CombinedGLM_Choice_3mm_Multimodel_Searchlight_%dmm_%dss_Level2/%s/%s_%s',...
                searchlight_radius,step_size,method_str,condition_str,modelNames{RDM});
            cur_folder = fullfile(scriptPath, secondLevelFolder);    


            % navigate to the relevant analysis directory 
            groupImFile = fullfile(cur_folder, 'SPM.mat');
            matlabbatch{1}.spm.stats.con.spmmat = {groupImFile};
            matlabbatch{1}.spm.stats.con.consess{1}.tcon.name = 'specified';
            matlabbatch{1}.spm.stats.con.consess{1}.tcon.weights = 1;
            matlabbatch{1}.spm.stats.con.consess{1}.tcon.sessrep = 'none';
            matlabbatch{1}.spm.stats.con.consess{2}.tcon.name = 'reverse';
            matlabbatch{1}.spm.stats.con.consess{2}.tcon.weights = -1;
            matlabbatch{1}.spm.stats.con.consess{2}.tcon.sessrep = 'none';
            matlabbatch{1}.spm.stats.con.delete = 0;
            spm('defaults','FMRI');
            spm_jobman('run', matlabbatch);   
        
        end 
    end 
 end 

cd(scriptPath);
clear all


