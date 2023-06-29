%% Extract the individual beta coefficents values for each ROI 

scriptPath = pwd; 
conPath = fullfile(scriptPath, 'Simple_DM_Timings_Level_2','Set1_Basic');
roiPath = fullfile(scriptPath,'ROI');

% load in roi list
roi_list = readtable('roi_list.xlsx');

% load in contrast list
contrast_list = readtable('contrasts_beta_extraction.xlsx');

%% Loop over each contrast and pull out betas for each ROI

% initialize an empty table to store data 
beta_table = {}; 

for con = 1:2
    
group = contrast_list.group(con); 
con_name = contrast_list.contrast(con);

    contrast = fullfile(conPath,con_name,'SPM.mat'); 
    table = []; 
    % loop over each of our ROIs
    for roi_number = 1:height(roi_list)

        % get name of roi 
        region = roi_list.roi_region(roi_number);
        table = create_roi_cols(contrast,region,con_name,scriptPath,roiPath); 
        beta_table = [beta_table; table]; 
    
    end 

end 


%% Save beta table         

filename = 'delta_beta_table.xlsx';                                   
writetable(beta_table,filename,'WriteVariableNames',true);   



