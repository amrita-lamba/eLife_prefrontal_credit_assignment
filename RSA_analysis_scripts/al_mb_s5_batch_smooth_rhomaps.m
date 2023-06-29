%% Code for smoothing whole brain maps with a 6mm smoothing kernal

subInfo = readtable('Adult_Run_List.csv');
scriptPath = pwd; 

% define searchlight radius
searchlight_radius = 9; 
step_size = 1;

% define the prediction RDM names here
condition = 1; 
subNum = 4002; % just grab from a random participant to get field names 
[pred_RDM,conditionVec] = construct_pred_RDMs_feedback(subNum,condition,scriptPath); % make sure you're grabbing the right matrices 
modelNames = fieldnames(pred_RDM);

%% Smooth brain maps 

% find scans to smooth
for sub = 1:30 
    
     subNum = subInfo.subNum(sub);
    
    if subInfo.exclusionCode(sub) == 0 
    
        for condition = 1:2

            % intialize empty scan vector
            scans = {};
            i = 0; 
            
             for RDM = 1:length(modelNames)

                % increment counter 
                i = i+1;   

                % navigate to subject's directory
                niftiPath = sprintf('%s/CombinedGLM_Feedback_3mm_Multimodel_Searchlight_%dmm_%dss/%d_searchlight_condition_%d',...
                            scriptPath,searchlight_radius,step_size,subNum,condition);

                conFile = dir(fullfile(niftiPath, '*.nii'));
                scans = [scans; ...
                reshape(strcat(niftiPath, '/', {conFile.name}, ',1'), [], 1)];      
                
             end
            
                matlabbatch = {};
                matlabbatch{1}.spm.spatial.smooth.data = scans; 
                matlabbatch{1}.spm.spatial.smooth.fwhm = [6 6 6]; % smoothing kernal
                matlabbatch{1}.spm.spatial.smooth.dtype = 0;
                matlabbatch{1}.spm.spatial.smooth.im = 0;
                matlabbatch{1}.spm.spatial.smooth.prefix = 's_';

                % run the batch
                spm_jobman('run', matlabbatch);
                clear matlabbatch
        end
    end
end



 



