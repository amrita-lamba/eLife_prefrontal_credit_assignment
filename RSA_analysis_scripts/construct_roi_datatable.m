%% read in data

%read subject info table with counter-balance and exclusion data 
subInfo = readtable('Adult_Run_List.csv');
scriptPath = pwd; 

% SPECIFY MODEL HERE & CONDITION HERE
RDM = 'identity_RDM';
RDM_roi = 1; % 1 = identity, 2 = qvalue
filename = 'feedback_identity_rsa_roi.xlsx'; 

% define user options here 
userOptions = defineUserOptions_socialRisk(scriptPath);

% grab the model names from any subject 
subNum = 4002;
condition = 1; 
[pred_RDM,conditionVec] = construct_pred_RDMs_choice(subNum,condition,scriptPath);
modelNames = fieldnames(pred_RDM);

% define search light parameters 
searchlightRad_mm = userOptions.searchlightRadius;
ss = 1; 

% initialize data columns
subj_col = [];
condition_col = []; 
roi_col = [];
model_col = []; 
rho_col = []; 
t_col = []; 

for sub = 1:30
    
    for condition = 1:2
    
        % subjects that will be included in final analysis 
        if subInfo.exclusionCode(sub) == 0 

        % extract out the relevant onset info for both runs
        subNum = subInfo.subNum(sub); 

        % navigate to subject directory
        roi_folder = (sprintf('%s/CombinedGLM_Feedback_3mm_Multimodel_Searchlight_%dmm_%dss/%d_searchlight_condition_%d/',...
        scriptPath,searchlightRad_mm, ss,subNum,condition));    
        
        % load in ROI data
        roi_mat = sprintf('%d_condition_%d_roi_coeffs_%s.mat',subNum,condition,RDM);
        load([roi_folder,roi_mat]);
        roi_names = fieldnames(roi_coeffs);
        roi_names_short = {};
        
        % remove binary 
        for roi = 1:length(roi_names)
            
            long_name = roi_names{roi};
            %short_name = extractBefore(string(long_name),"_binary");
            %roi_names_short{roi} = char(short_name);
            roi_names_short{roi} = char(long_name);
               
        end 

            subj_vals = repelem(subNum,length(roi_names),1)';
            condition_vals = repelem(condition,length(roi_names),1)';
            model_vals = repelem(modelNames(RDM_roi),length(roi_names),1);
            
            for roi = 1:length(roi_names)
                
                % check whethehr subject has data in ROI            
                if ~isempty(roi_coeffs.(roi_names{roi})) % we have data
                
                rho_val = roi_coeffs.(roi_names{roi}).b_est;
                rho_col = [rho_col, rho_val]; 
                
                t_val = roi_coeffs.(roi_names{roi}).t_stat;
                t_col = [t_col, t_val]; 
                
                else % no data
                    
                rho_val = NaN; 
                rho_col = [rho_col, rho_val]; 
                
                t_val = NaN; 
                t_col = [t_col, t_val]; 
                
                end 

            end 

            subj_col = [subj_col, subj_vals];
            condition_col = [condition_col, condition_vals];
            model_col = [model_col; model_vals];
            roi_col = [roi_col; roi_names_short'];

        elseif subInfo.exclusionCode(sub) == 1
            fprintf('Invalid Subject'); 
        end 

    end 
    
end 

% construct the data table 
table_array = subj_col'; 
data_table = array2table(table_array);
data_table.condition = condition_col';
data_table.model = model_col; 
data_table.roi = roi_col; 
data_table.rho = rho_col'; 
data_table.t = t_col'; 

% rename column 1 
data_table.Properties.VariableNames = {'subjNum','condition','model','roi','b_est','t_stat'}; 

% write the data table
writetable(data_table,filename,'WriteVariableNames',true);   



