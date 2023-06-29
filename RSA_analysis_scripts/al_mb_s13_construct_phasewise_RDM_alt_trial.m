%% Run the cross-timepoint RSA (choice + feedback) 

subInfo = readtable('Adult_Run_List.csv');
data = readtable('adult_behavData_2ca2lr.csv');
scriptPath = pwd; 

roiPath = [scriptPath,'/ROI_combinedGLM','/conjunction_identity']; % change this for each model
roiList = dir(fullfile(roiPath,'*.nii'));

%% load in feedback & choice betas

% IMPORTANT: specify your output directory here
fb_folder = [pwd,'/Trialwise_Voxel_Betas/Feedback_combined_GLM'];
choice_folder = [pwd,'/Trialwise_Voxel_Betas/Choice_combined_GLM'];
output_folder = [pwd,'/Trialwise_Voxel_Betas/Phasewise_RDMs_Alt'];

%% loop over subjects

% construct hypothesis matrices
id_mat = (eye(4)*-1)+1;
ones_mat = ones(4,4);
zeros_mat = zeros(4,4); 

h1 = [id_mat,ones_mat; ones_mat,ones_mat];
h2 = [ones_mat,ones_mat; ones_mat,id_mat];
h3 = [ones_mat,id_mat; id_mat,ones_mat];

% create logic vector of what to include in each matrix
h1_logic = [ones_mat,zeros_mat; zeros_mat, zeros_mat];
h2_logic = [zeros_mat,zeros_mat; zeros_mat, ones_mat];
h3_logic = [zeros_mat,ones_mat; ones_mat, zeros_mat];

% turn into vector
h1_vec = zscore(h1(:));
h2_vec = zscore(h2(:));
h3_vec = zscore(h3(:));

% turn into vector
h1_logic_vec = h1_logic(:);
h2_logic_vec = h2_logic(:);
h3_logic_vec = h3_logic(:);

% select only phases we care about 
h1_idx = find(h1_logic_vec ~= 0);
h2_idx = find(h2_logic_vec ~= 0);
h3_idx = find(h3_logic_vec ~= 0);

for subIdx = 1:30
    
    roi_mat = struct; 

    % subjects that will be included in final analysis 
    if subInfo.exclusionCode(subIdx) == 0 
        
        subNum = subInfo.subNum(subIdx); 
        
        % loop over task
        for condition = 1:2
        
            % pull voxel activations from each ROI for each phase
            fb_betas = load(sprintf('%s/%d_condition_%d_voxel_actvs',...
            fb_folder,subNum,condition));

            choice_betas = load(sprintf('%s/%d_condition_%d_voxel_actvs',...
            choice_folder,subNum,condition));
        
            % construct & save the matrix
            sub_dat = data(data.subjNum == subNum & data.condition == condition,:); 
            
            % loop over each ROI
            for roi = 1:length(roiList)
                
                file_str = roiList(roi).name;  
                roi_name = extractBefore(string(file_str),".nii");

                % initialize roi stats component of structure
                roi_mat.(roi_name) = [];    

                % initialize matrices
                dim = size(fb_betas.roi_coeffs.(roi_name).trial_1.vox_betas);
                even_fb_mat = zeros(dim(2),60);
                even_choice_mat = zeros(dim(2),60);
                odd_fb_mat = zeros(dim(2),60);
                odd_choice_mat = zeros(dim(2),60);
                
                even_fb_player_mat = zeros(60,1);
                even_choice_player_mat = zeros(60,1);
                odd_fb_player_mat = zeros(60,1);
                odd_choice_player_mat = zeros(60,1);
               
                for i = 1:height(sub_dat)

                    binned_trial = sub_dat.binned_trial_num(i);
                    pNum = sub_dat.playerNum(i);
                    trial_string = sprintf('trial_%d',binned_trial);
                    trial_div = rem(binned_trial,2);
                    
                    if trial_div == 0  
                        
                        even_choice_mat(:,i) = choice_betas.roi_coeffs.(roi_name).(trial_string).vox_betas(pNum,:);
                        even_choice_player_mat(i) = pNum;
                        even_fb_mat(:,i) = fb_betas.roi_coeffs.(roi_name).(trial_string).vox_betas(pNum,:);
                        even_fb_player_mat(i) = pNum;
                        
                    else % odd trials
                        
                        odd_choice_mat(:,i) = choice_betas.roi_coeffs.(roi_name).(trial_string).vox_betas(pNum,:);
                        odd_choice_player_mat(i) = pNum;
                        odd_fb_mat(:,i) = fb_betas.roi_coeffs.(roi_name).(trial_string).vox_betas(pNum,:);
                        odd_fb_player_mat(i) = pNum;
                    end 
                   
                end 
                
                % store matrices
                roi_mat.(roi_name).even_fb_mat = even_fb_mat;     
                roi_mat.(roi_name).even_choice_mat = even_choice_mat;
                roi_mat.(roi_name).odd_fb_mat = odd_fb_mat;     
                roi_mat.(roi_name).odd_choice_mat = odd_choice_mat;
                
                % calculate correlation distances & compute RDMs
                % first, compute mean activations for each player type
                even_mean_choice_actvs = zeros(dim(2),4);
                even_mean_fb_actvs = zeros(dim(2),4);
                odd_mean_choice_actvs = zeros(dim(2),4);
                odd_mean_fb_actvs = zeros(dim(2),4);
                
                for p = 1:4
                    
                    even_mean_choice_actvs(:,p)=mean(even_choice_mat(:,even_choice_player_mat==p), 2);
                    even_mean_fb_actvs(:,p)=mean(even_fb_mat(:,even_fb_player_mat==p), 2);
                    odd_mean_choice_actvs(:,p)=mean(odd_choice_mat(:,odd_choice_player_mat==p), 2);
                    odd_mean_fb_actvs(:,p)=mean(odd_fb_mat(:,odd_fb_player_mat==p), 2);

                end

                even_all_actvs = [even_mean_choice_actvs, even_mean_fb_actvs];
                odd_all_actvs = [odd_mean_choice_actvs, odd_mean_fb_actvs];
                
                % compuate correlation distance 
                clear corrDist
                corrDist=1-corr(even_all_actvs,odd_all_actvs); % this is matrix produce
                z_corr_d = zscore(corrDist);
%                 scriptPath(corrDist)
%                 colorbar()
                
                %compute euclidean distance
                clear eucDist
                for i = 1:size(even_all_actvs, 2)
                    for j = 1:size(even_all_actvs, 2)
                       eucDist(i,j)= sqrt(sum((even_all_actvs(:,i)-odd_all_actvs(:,j)).^2));
                    end
                end
%                 scriptPath(eucDist)
%                 colorbar()
                               
                % save all stats
                roi_mat.(roi_name).even_mean_choice_actvs = even_mean_choice_actvs;
                roi_mat.(roi_name).odd_mean_choice_actvs = odd_mean_choice_actvs;
                roi_mat.(roi_name).even_mean_fb_actvs = even_mean_fb_actvs; 
                roi_mat.(roi_name).odd_mean_fb_actvs = odd_mean_fb_actvs; 
                roi_mat.(roi_name).even_all_actvs = even_all_actvs; 
                roi_mat.(roi_name).odd_all_actvs = odd_all_actvs; 
                roi_mat.(roi_name).corrDist = corrDist; 
                
                % reshape matrices
                corrDist_vec = z_corr_d(:);
                eucDist_vec = eucDist(:);
                
                % control for all other regressors (more conservative)
%                 data_table = [corrDist_vec,h1_vec,h2_vec,h3_vec];
%                 data_table = array2table(data_table);
%                 data_table.Properties.VariableNames = {'corrDist','h1','h2','h3'}; 
%                 lm1 = fitlm(data_table,'corrDist~h1');
%                 lm2 = fitlm(data_table,'corrDist~h2');
%                 lm3 = fitlm(data_table,'corrDist~h3');

                % save 
                roi_mat.(roi_name).choice_identity = corr(corrDist_vec(h1_idx),h1_vec(h1_idx));
                roi_mat.(roi_name).fb_identity = corr(corrDist_vec(h2_idx),h2_vec(h2_idx)); 
                roi_mat.(roi_name).cross_identity = corr(corrDist_vec(h3_idx),h3_vec(h3_idx)); 

            end 
            
            % save coefficients
            mat_file_name = (sprintf('%s/%d_condition_%d_phasewise_rdms',...
            output_folder,subNum,condition));
            save(mat_file_name, 'roi_mat');
            
        end     
    end     
end

