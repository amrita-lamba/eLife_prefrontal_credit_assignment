%% Combine the TG & SM beta maps 

%% Note: Make sure the single subject/condition-wise beta images are smoothed first 

%% combine images 

subInfo = readtable('Adult_Run_List.csv');
scriptPath = pwd; 
combinedPath = [pwd,'/CombinedGLM_Choice_3mm_Multimodel_Searchlight_9mm_1ss'];
rsaPath = [pwd,'/CombinedGLM_Choice_3mm_Multimodel_Searchlight_9mm_1ss_Level2'];

method_type = {'b_est','corr_dist','corr_dist_ztrans','t_stat'};

% define searchlight radius
searchlight_radius = 9; 
step_size = 1;

% define the prediction RDM names here
condition = 1; 
subNum = 4002; % just grab from a random participant to get field names 
[pred_RDM,conditionVec] = construct_pred_RDMs_choice(subNum,condition,scriptPath);  % make sure you're grabbing the right matrices 
modelNames = fieldnames(pred_RDM);

for sub = 1:30 
    
    subNum = subInfo.subNum(sub);
    
    if subInfo.exclusionCode(sub) == 0 
        
        for RDM = 1:length(modelNames)
                
                for m = 1:length(method_type)
                    
                method_str = char(method_type(m));
                               
                % find the relevant .nii images for each RDM & each
                % condition
                scans = {}; 
                cnt = 0; 
                    for condition = 1:2
                        
                        % pull the appropriate level 2 tmap (binarized)
                        if condition == 1
                        condition_str = 'TG';                
                        else
                        condition_str = 'SM';
                        end 
                        
                        cnt = cnt + 1; 
                        sub_beta_folder = sprintf('/%d_searchlight_condition_%d',subNum,condition);
                        beta_path = [combinedPath,sub_beta_folder];
                        nii_file = sprintf('/s_sub-%d_%s_%s.nii',subNum,modelNames{RDM},method_str); % pull smoothed image
                        scans{cnt} = [beta_path,nii_file];
                        
                    end 
                    
                % create subject output directory for combined beta maps 
                subject_combined_folder = sprintf('/%d_searchlight_combined',subNum);
                outputPath = [combinedPath,subject_combined_folder];
                
                %make the directory if it doesn't exist
                if ~exist(outputPath)
                     mkdir(outputPath)
                end
                    
                % combine scans 
                matlabbatch = {}; 
                matlabbatch{1}.spm.util.imcalc.input = scans'; 
                matlabbatch{1}.spm.util.imcalc.output = sprintf('s_sub-%d_%s_%s_combined',subNum,modelNames{RDM},method_str);
                matlabbatch{1}.spm.util.imcalc.outdir = {outputPath};
                matlabbatch{1}.spm.util.imcalc.expression = 'i1+i2'; % use summed map (tg + sm)
                matlabbatch{1}.spm.util.imcalc.var = struct('name', {}, 'value', {});
                matlabbatch{1}.spm.util.imcalc.options.dmtx = 0;
                matlabbatch{1}.spm.util.imcalc.options.mask = 0;
                matlabbatch{1}.spm.util.imcalc.options.interp = 1;
                matlabbatch{1}.spm.util.imcalc.options.dtype = 4;
                
                % run the batch
                spm_jobman('run', matlabbatch); 
                clear matlabbatch
                
               end 
        end 
    end 
end 

cd(scriptPath);
clear all

%% Run the second level analysis on combined beta maps 

subInfo = readtable('Adult_Run_List.csv');
scriptPath = pwd; 
combinedPath = [pwd,'/CombinedGLM_Choice_3mm_Multimodel_Searchlight_9mm_1ss'];
rsaPath = [pwd,'/CombinedGLM_Choice_3mm_Multimodel_Searchlight_9mm_1ss_Level2'];

method_type = {'b_est','corr_dist','corr_dist_ztrans','t_stat'};

% define searchlight radius
searchlight_radius = 9; 
step_size = 1;

% define the prediction RDM names here
condition = 1; 
subNum = 4002; % just grab from a random participant to get field names 
[pred_RDM,conditionVec] = construct_pred_RDMs_choice(subNum,condition,scriptPath);  % make sure you're grabbing the right matrices 
modelNames = fieldnames(pred_RDM);

for RDM = 1:length(modelNames)

        for m = 1:length(method_type)

        method_str = char(method_type(m));
        
        % navigate to the relevant analysis directory 
        cur_folder = sprintf('%s/%s/%s',secondLevelFolder, method_str,modelNames{RDM});
        
        if ~exist(cur_folder)
            mkdir(cur_folder)
        end

        % intialize empty scan vector
        scans = {};
        matlabbatch = {};
        matlabbatch{1}.spm.stats.factorial_design.dir = {cur_folder};

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
                          outputPath,subNum,modelNames{RDM},method_str); % smoothed/unsmoothed

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
        matlabbatch{1}.spm.stats.factorial_design.masking.em = {'/gpfs/data/ofeldman/alamba1/adult_social_risk/Final_pipeline_scripts/Final_RSA_analysis_scripts/Perm_Masks/conj_identity_tstat_perm_05_binary.nii'};
        matlabbatch{1}.spm.stats.factorial_design.globalc.g_omit = 1;
        matlabbatch{1}.spm.stats.factorial_design.globalm.gmsca.gmsca_no = 1;
        matlabbatch{1}.spm.stats.factorial_design.globalm.glonorm = 1;
        spm_jobman('run', matlabbatch); 

        end 
      
end 

cd(scriptPath);
clear all

%% Estimate the combined beta maps 

scriptPath = pwd; 
combinedPath = [pwd,'/CombinedGLM_Choice_3mm_Multimodel_Searchlight_9mm_1ss'];
secondLevelFolder = [pwd,'/Conj_CombinedGLM_Choice_3mm_Multimodel_Searchlight_9mm_1ss_Level2'];

method_type = {'b_est','corr_dist','corr_dist_ztrans','t_stat'};

% define searchlight radius
searchlight_radius = 9; 
step_size = 1;

% define the prediction RDM names here
condition = 1; 
subNum = 4002; % just grab from a random participant to get field names 
[pred_RDM,conditionVec] = construct_pred_RDMs_choice(subNum,condition,scriptPath);
modelNames = fieldnames(pred_RDM);

for RDM = 1:length(modelNames)

    % loop over each distance method
    for m = 1:length(method_type)
        
    method_str = char(method_type(m));

    % navigate to the relevant analysis directory 
    cur_folder = sprintf('%s/%s/%s',secondLevelFolder, method_str,modelNames{RDM}); 

    % navigate to the relevant analysis directory 
    groupImFile = fullfile(cur_folder, 'SPM.mat');
    matlabbatch = {};
    matlabbatch{1}.spm.stats.fmri_est.spmmat = {groupImFile};
    matlabbatch{1}.spm.stats.fmri_est.write_residuals = 0;
    matlabbatch{1}.spm.stats.fmri_est.method.Classical = 1;
    spm_jobman('run', matlabbatch); 

    end 
end 

cd(scriptPath);
clear all 

%% T-test against 0

scriptPath = pwd; 
combinedPath = [pwd,'/CombinedGLM_Choice_3mm_Multimodel_Searchlight_9mm_1ss'];
rsaPath = [pwd,'/CombinedGLM_Choice_3mm_Multimodel_Searchlight_9mm_1ss_Level2'];

method_type = {'b_est','corr_dist','corr_dist_ztrans','t_stat'};

% define searchlight radius
searchlight_radius = 9; 
step_size = 1;

% define the prediction RDM names here
condition = 1; 
subNum = 4002; % just grab from a random participant to get field names 
[pred_RDM,conditionVec] = construct_pred_RDMs_choice(subNum,condition,scriptPath);
modelNames = fieldnames(pred_RDM);
    
for RDM = 1:length(modelNames)

    % loop over each distance method
    for m = 1:length(method_type)
        
    method_str = char(method_type(m));

    % navigate to the relevant analysis directory 
    cur_folder = sprintf('%s/%s/%s',secondLevelFolder, method_str,modelNames{RDM});

    % navigate to the relevant analysis directory 
    groupImFile = fullfile(cur_folder, 'SPM.mat');
    matlabbatch = {};
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

cd(scriptPath);
clear all
