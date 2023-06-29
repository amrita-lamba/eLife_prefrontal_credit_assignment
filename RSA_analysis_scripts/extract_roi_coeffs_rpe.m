function extract_roi_coeffs_rpe(subNum,scriptPath,pred_RDM,userOptions,condition,conditionVec,roiPath,roiList,output_folder,phase)

% initalize a structure to save rho and beta coeffsficients for each ROI
roi_coeffs = struct; 

% add the RSA toolbox to path
toolboxRoot = [scriptPath,'/rsatoolbox-develop/'];
addpath(genpath(toolboxRoot));
% addpath(genpath('/gpfs/data/ofeldman/alamba1/spm12/toolbox/marsbar'));

% GLM folder name
glmpath = [scriptPath,'/Full_DM_Timings_Combined/'];

% specify brain dims
xdim = 79;
ydim = 95;
zdim = 79;

% noise normalization options - overall will carry out noise normalization
% over whole dataset rather than run-wise
opt.normmode='overall';

% create model RDMs
rawModels = pred_RDM; 
Models = rsa.constructModelRDMs(rawModels, userOptions);
M0 = [];
M0_notzscored = [];

%% Searchlight Analysis %% 
    
% get the spm file for this subject
timePath = [scriptPath,'/Full_DM_Timings_Combined/'];
sub = char(sprintf("%d_timings", subNum));
spm_file = fullfile(timePath, sub, 'SPM.mat');
load(spm_file);

% get subject brain mask
maskroiPath = [scriptPath,'/Subject_Masks/'];
roi = load(fullfile(maskroiPath, sprintf('wholebrain_%d.mat',subNum)));

% make marsbar design object from the SPM.mat containing bold info
design_obj = mardo(spm_file);
data_obj = get_marsy(roi.roi,design_obj,'mean');

% get summary time course(s) 
y = region_data(data_obj);
Y=y{1,1}; % scans/volumes x voxels 
xyz_coords = xyz(data_obj,1,'vox'); 

% read in whole brain mask volumes 
maskPath = [timePath, sub];
mask_header = spm_vol(fullfile(maskPath, 'mask.nii'));
mask_vol = spm_read_vols(mask_header);

% this part is coverting the time course Y into a 4D object
Y2 = zeros(size(Y,1),xdim,ydim,zdim);
for vi = 1:length(xyz_coords)
    x = xyz_coords(1,vi);
    y = xyz_coords(2,vi);
    z = xyz_coords(3,vi);
    Y2(:,x,y,z) = Y(:,vi);
        if mod(vi,1000) == 0
        end
end
Y2(Y2==0) = nan;

% calculate correlation coeffsficients for each model within that ROI
for k = 1:length(Models)

    % pull out each ROI
    for rsa_roi = 1:length(roiList)
    
        % get the roi name
        file_str = roiList(rsa_roi).name;  
        roi_name = extractBefore(string(file_str),".nii");

        % initialize roi stats component of structure
        roi_coeffs.(roi_name) = [];
        file_name = roiList(rsa_roi).name;
        masktemp_header = spm_vol(fullfile(roiPath,file_name));
        masktemp = spm_read_vols(masktemp_header);
        masktemp = masktemp.*mask_vol;

        % pull out betas within the ROI     
        Ytemp = Y2(:,logical(masktemp));
        bid = find(~isnan(sum(Ytemp)));
        dim = size(Ytemp);
        
        % does participant have brain data in this region?
        if ~isempty(Ytemp) && dim(2) > 2
            
            [betas] = noiseNormalizeBeta_2015(Ytemp(:,bid),SPM,opt);

            % calculate distances using different methods
            rpe_names = fieldnames(conditionVec);
            rpe_string = char(rpe_names(k));
            corr_d = calc_neural_corr_dist(betas,conditionVec.(rpe_string)(:,1));
            euc_d = calc_neural_euc_dist(betas,conditionVec.(rpe_string)(:,1));

            % get vectors of model predicted dis-similarities and z-score them
            M0(1:length(corr_d),k) = zscore(Models(k).RDM(find(tril(ones(size(Models(k).RDM)),-1))));
            M0_notzscored(1:length(corr_d),k) = Models(k).RDM(find(tril(ones(size(Models(k).RDM)),-1))); 

            % normalize 
            z_corr_d = zscore(corr_d);

            % get the fieldames & store distance values 
            modelnames = fieldnames(rawModels);
            roi_coeffs.(roi_name).corr_d = corr_d; 
            roi_coeffs.(roi_name).z_corr_d = z_corr_d; 

            % calculate similarity between model and neural data 
            corr_rho = corr(z_corr_d,M0(1:length(corr_d),k),'type','Pearson');
       
            % z transform the rho coeffsficents
            corr_rho_ztrans =.5.*log((1+corr_rho)./(1-corr_rho));

            % populate matrix
            V = euc_d-M0(1:length(euc_d));
            D = sqrt(V * V');
            euc = D;
           
            % store data
            roi_coeffs.(roi_name).rho = corr_rho;
            roi_coeffs.(roi_name).fisher_z_rho = corr_rho_ztrans;
            roi_coeffs.(roi_name).euc = euc;
           
        end
        
    end 
    
    % save .mat file
    mat_file_name = (sprintf('%s/%d_searchlight_condition_%d/%d_condition_%d_roi_coeffs_2ca2lr_zDelta_%s_%s',...
    output_folder,subNum,condition,subNum,condition,modelnames{k},phase));
    save(mat_file_name, 'roi_coeffs');
    

end 

end 

