%% Run the second level analyses for group level analyses

%% Create the group level .SPM file for basic contrasts
subInfo = readtable('Adult_Run_List.csv');
con_list = readtable('contrast_masterList.csv');
scriptPath = pwd;
conPath = fullfile(scriptPath, 'Simple_DM_Timings_Level_2','Set1_Basic');

for contrast = 1:10

    A1 = con_list.conNumber(contrast);
    A2 = con_list.conName(contrast);

    A2 = char(A2);

    contrast_number = contrast; 

    con_folder = '%d_%s'; 
    con_folder_name = sprintf(con_folder,A1,A2);

    % navigate to the relevant analysis directory 
    cur_folder = fullfile(conPath, con_folder_name);    
    matlabbatch{1}.spm.stats.factorial_design.dir = {cur_folder};

    scans = {}; 
    timePath = fullfile(scriptPath,'Simple_DM_Timings');

    %find that subject's .con file
    i = 0; 
    for sub = 1:height(subInfo)
    
        subNum = subInfo.subNum(sub);
        A1 = subNum;

        if subInfo.exclusionCode(sub) == 0 
        i = i + 1;     
        str_subj = "%d_timings";
        subTimings = sprintf(str_subj, A1);
        subTimings = char(subTimings);
        cd(fullfile(timePath, subTimings));    
        conFile = dir(fullfile(timePath, subTimings, 'con_*.nii'));
        conFile = struct2cell(conFile);
        conFile = conFile(1,contrast);
        conFile = cell2mat(conFile);

        conFilePath = fullfile(timePath, subTimings,conFile);
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
    matlabbatch{1}.spm.stats.factorial_design.masking.em = {''};
    matlabbatch{1}.spm.stats.factorial_design.globalc.g_omit = 1;
    matlabbatch{1}.spm.stats.factorial_design.globalm.gmsca.gmsca_no = 1;
    matlabbatch{1}.spm.stats.factorial_design.globalm.glonorm = 1;

    spm_jobman('run', matlabbatch); 
      
end 

cd(scriptPath);
clear all

%% 2) Estimate the model  

% enter name of .csv file with contrast to run in here 
con_list = readtable('contrast_masterList.csv');
scriptPath = pwd;
conPath = fullfile(scriptPath, 'Simple_DM_Timings_Level_2','Set1_Basic');

for contrast = 1:10

    A1 = con_list.conNumber(contrast);
    A2 = con_list.conName(contrast); 
    A2 = char(A2);

    con_folder = '%d_%s'; 
    con_folder_name = sprintf(con_folder,A1,A2);

    % navigate to the relevant analysis directory 
    groupImFile = fullfile(conPath, con_folder_name, 'SPM.mat');

    matlabbatch{1}.spm.stats.fmri_est.spmmat = {groupImFile};
    matlabbatch{1}.spm.stats.fmri_est.write_residuals = 0;
    matlabbatch{1}.spm.stats.fmri_est.method.Classical = 1;

    spm_jobman('run', matlabbatch);     
end 

cd(scriptPath);
clear all 

%% 3) Run the relevant one-sample t-tests 

% enter name of .csv file with contrast to run in here 
con_list = readtable('contrast_masterList.csv');
scriptPath = pwd;
conPath = fullfile(scriptPath, 'Simple_DM_Timings_Level_2','Set1_Basic');

for contrast = 1:10
    
    A1 = con_list.conNumber(contrast);
    A2 = con_list.conName(contrast);
    A2 = char(A2);

    con_folder = '%d_%s'; 
    con_folder_name = sprintf(con_folder,A1,A2);

    % navigate to the relevant analysis directory 
    groupImFile = fullfile(conPath, con_folder_name, 'SPM.mat');    

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

cd(scriptPath);
clear all

