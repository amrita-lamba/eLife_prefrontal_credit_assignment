%% this script will create an ROI mask for each subject
function create_individual_masks(subNum,scriptPath) 

% define path for subject's individual brain mask
timePath = [scriptPath,'/Full_DM_Timings_Combined/'];
sub = char(sprintf("%d_timings", subNum));
imgname = fullfile(timePath, sub, 'mask.nii');

% convert nifti to roi object
o = maroi_image(struct('vol', spm_vol(imgname), 'binarize',0,...
'func', 'img'));
filename = sprintf("wholebrain_%d", char(subNum));

% save file in mask directory 
cd('Subject_Masks')
saveroi(o, filename);
cd(scriptPath)
