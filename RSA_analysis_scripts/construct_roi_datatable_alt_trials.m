%% write cross-correlations from even-odd analysis to a datatable 

scriptPath = pwd; 
subInfo = readtable('Adult_Run_List.csv');
roiPath = [scriptPath,'/ROI_combinedGLM','/conjunction_identity']; % change this for each model
roiList = dir(fullfile(roiPath,'*.nii'));
output_folder = [pwd,'/Trialwise_Voxel_Betas/Phasewise_RDMs_Alt'];

%% loop over subjects

% initialize a structure to store data
roi_coef_TG = [];
roi_coef_SM = [];

% loop over condition 
for condition = 1:2

    % loop over ROIs
    for roi = 1:length(roiList)

        matrix = zeros(28,1);

        % initialize element counter
        idx = 0; 
        for subIdx = 1:30

          % subjects that will be included in final analysis 
           if subInfo.exclusionCode(subIdx) == 0 

            subNum = subInfo.subNum(subIdx); 
            idx = idx + 1; 

            % load in data 
            mat_file_name = (sprintf('%s/%d_condition_%d_phasewise_rdms',...
            output_folder,subNum,condition));
            load(mat_file_name);

            file_str = roiList(roi).name;  
            roi_name = extractBefore(string(file_str),".nii");

            matrix(idx,1) = roi_mat.(roi_name).choice_identity;
            matrix(idx,2) = roi_mat.(roi_name).fb_identity;
            matrix(idx,3) = roi_mat.(roi_name).cross_identity;
            
            % update TG and SM matrices 
            if condition == 1
                roi_coef_TG(:,:,roi) = matrix; 
            else
                roi_coef_SM(:,:,roi) = matrix;     
            end 

            end    
        end     
    end
end 

%% write the data table 

% SPECIFY MODEL HERE & CONDITION HERE
filename = 'alt_conj_rho.xlsx'; 

% initialize data columns
subj_col = [];
condition_col = []; 
h_col = {}; 
roi_col = {};
rho_col = []; 

% define predictors 
h_list = {'choice','feedback','crossphase'}';

% get subject list 
sub_list = table2array(subInfo(subInfo.exclusionCode == 0,2));
    
for condition = 1:2

    % loop over rois 
    for roi = 1:length(roiList)

        file_str = roiList(roi).name;  
        roi_name = extractBefore(string(file_str),".nii");

        if condition == 1
        roi_matrix = roi_coef_TG(:,:,roi);

        else 
        roi_matrix = roi_coef_SM(:,:,roi);

        end 

        % loop over predictors 
        for b = 1:3

            % h_vec
            h_vec = roi_matrix(:,b);
            
            % update all indicies 
            subj_col = [subj_col; sub_list];
            condition_col = [condition_col; repelem(condition,length(h_vec))']; 
            roi_col = [roi_col; repelem(roi_name,length(h_vec))'];
            h_col = [h_col;repelem(string(h_list{b}),length(h_vec))']; 
            rho_col = [rho_col; h_vec]; 

        end 
    end 
end 

% construct the data table 
table_array = subj_col; 
data_table = array2table(table_array);
data_table.condition = condition_col;
data_table.roi = roi_col;  
data_table.h = h_col;
data_table.rho = rho_col; 

% rename column 1 
data_table.Properties.VariableNames = {'subjNum','condition','roi','h','rho'}; 

% write the data table
writetable(data_table,filename,'WriteVariableNames',true);























