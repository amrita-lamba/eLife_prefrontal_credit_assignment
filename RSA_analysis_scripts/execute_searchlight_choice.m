function execute_searchlight_choice(subNum,scriptPath,pred_RDM,userOptions,condition,ss,conditionVec)

% Searchlight with model RDM for trial and reward in social vs nonsocial
% reward learning paradigm 
% code based on RSA script from Avinash Vaidya

%% Initialize  %%

% add the RSA toolbox to path
toolboxRoot = [scriptPath,'/rsatoolbox-develop/'];
addpath(genpath(toolboxRoot));
addpath(genpath('/gpfs/data/ofeldman/alamba1/spm12/toolbox/marsbar'));

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

% get vectors of model predicted dis-similarities and z-score them
for i = 1:length(Models)
    M0(:,i) = zscore(Models(i).RDM(find(tril(ones(size(Models(i).RDM)),-1))));
    M0_notzscored(:,i) = Models(i).RDM(find(tril(ones(size(Models(i).RDM)),-1))); 
end

%% Searchlight Analysis %% 
    
% get the spm file for this subject
timePath = [scriptPath,'/Full_DM_Timings_Combined/'];
sub = char(sprintf("%d_timings", subNum));
spm_file = fullfile(timePath, sub, 'SPM.mat');
load(spm_file);

% get group level brain mask
roiPath = [scriptPath,'/Subject_Masks/'];
roi = load(fullfile(roiPath, sprintf('wholebrain_%d.mat',subNum)));

% make marsbar design object from the SPM.mat containing bold info
design_obj = mardo(spm_file);
data_obj = get_marsy(roi.roi,design_obj,'mean');

% get summary time course(s) 
y = region_data(data_obj);
Y=y{1,1}; % scans/volumes x voxels 
xyz_coords = xyz(data_obj,1,'vox');

% read in mask volumes 
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

% create spherical multivariate searchlight - this bit of code is stolen from RSA
% Toolbox directly
searchlightRad_mm = userOptions.searchlightRadius;
voxSize_mm = userOptions.voxelSize;
rad_vox=searchlightRad_mm./voxSize_mm;
minMargin_vox=floor(rad_vox);
[x,y,z]=meshgrid(-minMargin_vox(1):minMargin_vox(1),-minMargin_vox(2):minMargin_vox(2),-minMargin_vox(3):minMargin_vox(3));
sphere=((x*voxSize_mm(1)).^2+(y*voxSize_mm(2)).^2+(z*voxSize_mm(3)).^2)<=(searchlightRad_mm^2);  % volume with sphere voxels marked 1 and the outside 0
sphereSize_vox=[size(sphere),ones(1,3-ndims(sphere))]; % enforce 3D (matlab stupidly autosqueezes trailing singleton dimensions to 2D, try: ndims(ones(1,1,1)). )
sphereLims = floor(sphereSize_vox/2);

% initialize empty mask for adding to sphere
emptymask = zeros(xdim,ydim,zdim);

% initialize empty brain volumes for each of model regressors
for b = 1:size(M0,2)
    statmaps(b).corr_dist = zeros(xdim,ydim,zdim);
    statmaps(b).corr_dist_ztrans = zeros(xdim,ydim,zdim);
    statmaps(b).b_est = zeros(xdim,ydim,zdim);
    statmaps(b).t_stat = zeros(xdim,ydim,zdim);
end

% initialize spotlight counts
si = 0;

% decide lower limit for volume to calculate distance on - if there are
% not enough voxels then don't try to calculate distances.
vollim = round(sum(sum(sum(sphere)))/5);

% now move the sphere around the brain
for xi = 1:ss:xdim
    for yi=1:ss:ydim
        for zi=1:ss:zdim
            % reset temp mask
            masktemp=emptymask;
            % get indecies for sphere
              xid = xi-sphereLims(1):xi+sphereLims(1);
              yid = yi-sphereLims(2):yi+sphereLims(2);
              zid = zi-sphereLims(3):zi+sphereLims(3);
              sphereTemp = sphere(xid > 0 & xid <= xdim,...
                                  yid > 0 & yid <= ydim,...
                                  zid > 0 & zid <= zdim);
              xid = xid(xid > 0 & xid <= xdim);
              yid = yid(yid > 0 & yid <= ydim);
              zid = zid(zid > 0 & zid <= zdim);
            % assign sphere to temp mask
            masktemp(xid,yid,zid)=sphereTemp;
            % mask the spot light with brain mask
            masktemp = masktemp.*mask_vol;
            % get brain data from mask
            Ytemp = Y2(:,logical(masktemp));
            % get indecies for voxels that are not outside brain
            bid = find(~isnan(sum(Ytemp)));
            % break  if no voxels have brain, or if there are too few
            % voxels
            if isempty(bid) || size(Ytemp,2) < vollim
                continue
            end

            % grab betas (n regressors from design matrix x n voxels in searchlight)
            [betas] = noiseNormalizeBeta_2015(Ytemp(:,bid),SPM,opt);
            
            % calculate distances using different methods
            corr_d = calc_neural_corr_dist(betas,conditionVec(:,1));
            
            % normalize 
            z_corr_d = zscore(corr_d);
       
            % store coefficients  
            for k = 1:length(statmaps)
                
                M0_table = array2table(M0);
                M0_table.Properties.VariableNames = fieldnames(rawModels);
                M0_table.z_corr_d = z_corr_d; 

                % calculate disimilarity between model and neural data 
                corr_rho = corr(z_corr_d,M0(:,k),'type','Pearson');
                
                % z transform the rho coefficents
                corr_rho_ztrans =.5.*log((1+corr_rho)./(1-corr_rho));
                
                % save stats
                statmaps(k).corr_dist(xi,yi,zi) = corr_rho;
                statmaps(k).corr_dist_ztrans(xi,yi,zi) = corr_rho_ztrans;
                
                if k == 1 
                    
                model_b = fitlm(M0_table,'z_corr_d~identity_RDM+trial_RDM');
                b_est = model_b.Coefficients.Estimate(2);
                t_stat = model_b.Coefficients.tStat(2);
                
                % save stats
                statmaps(k).b_est(xi,yi,zi) = b_est;
                statmaps(k).t_stat(xi,yi,zi) = t_stat;

                elseif k == 2
                    
                model_b = fitlm(M0_table,'z_corr_d~qval_RDM+trial_RDM');
                b_est = model_b.Coefficients.Estimate(2);
                t_stat = model_b.Coefficients.tStat(2);
                
                % save stats
                statmaps(k).b_est(xi,yi,zi) = b_est;
                statmaps(k).t_stat(xi,yi,zi) = t_stat;
                    
                elseif k == 3
                    
                model_b = fitlm(M0_table,'z_corr_d~trial_RDM');
                b_est = model_b.Coefficients.Estimate(2);
                t_stat = model_b.Coefficients.tStat(2);
                
                % save stats
                statmaps(k).b_est(xi,yi,zi) = b_est;
                statmaps(k).t_stat(xi,yi,zi) = t_stat;
                    
                end                   
            end
        end
    end
end

% save the brains for each subject
% load a volume (any volume) that is in the same  space as the searchlight analysis    
outputV = spm_vol(spm_vol(fullfile(timePath, sub, 'beta_0001.nii')));
outputV=rmfield(outputV,'pinfo');

% Specify where you want to save these data
outputDir = sprintf('%s/CombinedGLM_Choice_3mm_Multimodel_Searchlight_%dmm_%dss',scriptPath,searchlightRad_mm,ss);
rmpath(genpath('/gpfs/data/ofeldman/alamba1/spm12/toolbox/marsbar'));

% create folder to save subject data
subjFolder = sprintf('%s/%d_searchlight_condition_%d/',outputDir, subNum,...
condition);

if ~exist(subjFolder)
    mkdir(subjFolder)
end

mat_file_name = (sprintf('%s/CombinedGLM_Choice_3mm_Multimodel_Searchlight_%dmm_%dss/%d_searchlight_condition_%d/%d_condition_%d',...
scriptPath,searchlightRad_mm, ss,subNum,condition,subNum,condition));

% create the files
modelnames = fieldnames(rawModels);
for k = 1:length(statmaps)
    
    % save correlation distance 
    outputV.fname=sprintf('%s/%d_searchlight_condition_%d/sub-%d_%s_corr_dist.nii',outputDir, subNum,...
    condition,subNum, modelnames{k});
    volOut = spm_write_vol(outputV,statmaps(k).corr_dist);
    
    outputV.fname=sprintf('%s/%d_searchlight_condition_%d/sub-%d_%s_corr_dist_ztrans.nii',outputDir, subNum,...
    condition,subNum, modelnames{k});
    volOut = spm_write_vol(outputV,statmaps(k).corr_dist_ztrans);
   
    outputV.fname=sprintf('%s/%d_searchlight_condition_%d/sub-%d_%s_b_est.nii',outputDir, subNum,...
    condition,subNum, modelnames{k});
    volOut = spm_write_vol(outputV,statmaps(k).b_est);
    
    outputV.fname=sprintf('%s/%d_searchlight_condition_%d/sub-%d_%s_t_stat.nii',outputDir, subNum,...
    condition,subNum, modelnames{k});
    volOut = spm_write_vol(outputV,statmaps(k).t_stat);
    
end 

% Save .mat file
save(mat_file_name, 'statmaps');
end 
