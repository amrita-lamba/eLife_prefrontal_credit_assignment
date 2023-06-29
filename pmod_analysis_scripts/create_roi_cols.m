function [table] = create_roi_cols(contrast,region,con_name,scriptPath,roiPath)

sub_list = [4002,4003,4005,4006,4007,4008,4009,4010,4011,4012,4013,4014,4015,4016,4017,4018,4019,4020,...
                   4021,4022,4023,4024,4025,4026,4027,4028,4029,4030]; 
               
% create file name 
A1 = char(region); 
formatSpec = '%s_roi.mat';
file_name = sprintf(formatSpec,A1); 

% load in roi 
roi_file = dir(fullfile(roiPath, file_name));
cd(roiPath);
load(roi_file.name);
R = roi; 

% load contrast 
cd(scriptPath); 
design = mardo(contrast{1}); 

% get individual betas
[individual_b,pc] = extract_roi_betas(R, design); 

%% create a table columns 
region_col = repmat(region,length(sub_list),1); 
con_col = repmat(con_name,length(sub_list),1); 

% subject column 
subject_col = sub_list';

% beta values
betas = individual_b; 

% p values
pc_col = repelem(pc,length(sub_list),1); 

% create table
table = [subject_col,betas,pc_col]; 
table = array2table(table); 
table.brain_region = region_col; 
table.contrast = con_col; 
table.Properties.VariableNames = {'subject','beta','p_corrected','brain_region','contrast'};
 
