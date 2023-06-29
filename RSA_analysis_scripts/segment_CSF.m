function segment_CSF(subjnum,scriptPath)

% Default directories
subjPath = fullfile(scriptPath,'NIFTI/', int2str(subjnum));

mpragedir =  dir(fullfile(subjPath, '*mprage*'));
mpragedir = struct2cell(mpragedir);
mpragedir = mpragedir(1,:);

% **********************************************************************

cd(subjPath)

%% This section will run the spm segmentation function
anatPath = fullfile(mpragedir{1});
scanfnames = dir(fullfile(subjPath,anatPath, '4*.img'));
volImage = {}; 
volImage = [volImage; ...
    reshape(strcat(anatPath, '/', {scanfnames.name}, ',1'), [], 1)];  

% find anatomical 
matlabbatch = {};
matlabbatch{1}.spm.spatial.preproc.channel.vols = volImage; 
matlabbatch{1}.spm.spatial.preproc.channel.biasreg = 0.001;
matlabbatch{1}.spm.spatial.preproc.channel.biasfwhm = 60;
matlabbatch{1}.spm.spatial.preproc.channel.write = [0 0];
matlabbatch{1}.spm.spatial.preproc.tissue(1).tpm = {'/oscar/data/ofeldman/alamba1/spm12/TPM.nii,1'};
matlabbatch{1}.spm.spatial.preproc.tissue(1).ngaus = 1;
matlabbatch{1}.spm.spatial.preproc.tissue(1).native = [1 0];
matlabbatch{1}.spm.spatial.preproc.tissue(1).warped = [0 1];
matlabbatch{1}.spm.spatial.preproc.tissue(2).tpm = {'/oscar/data/ofeldman/alamba1/spm12/tpm/TPM.nii,2'};
matlabbatch{1}.spm.spatial.preproc.tissue(2).ngaus = 1;
matlabbatch{1}.spm.spatial.preproc.tissue(2).native = [1 0];
matlabbatch{1}.spm.spatial.preproc.tissue(2).warped = [0 1];
matlabbatch{1}.spm.spatial.preproc.tissue(3).tpm = {'/oscar/data/ofeldman/alamba1/spm12/tpm/TPM.nii,3'};
matlabbatch{1}.spm.spatial.preproc.tissue(3).ngaus = 2;
matlabbatch{1}.spm.spatial.preproc.tissue(3).native = [1 0];
matlabbatch{1}.spm.spatial.preproc.tissue(3).warped = [0 1];
matlabbatch{1}.spm.spatial.preproc.tissue(4).tpm = {'/oscar/data/ofeldman/alamba1/spm12/tpm/TPM.nii,4'};
matlabbatch{1}.spm.spatial.preproc.tissue(4).ngaus = 3;
matlabbatch{1}.spm.spatial.preproc.tissue(4).native = [1 0];
matlabbatch{1}.spm.spatial.preproc.tissue(4).warped = [0 1];
matlabbatch{1}.spm.spatial.preproc.tissue(5).tpm = {'/oscar/data/ofeldman/alamba1/spm12/tpm/TPM.nii,5'};
matlabbatch{1}.spm.spatial.preproc.tissue(5).ngaus = 4;
matlabbatch{1}.spm.spatial.preproc.tissue(5).native = [1 0];
matlabbatch{1}.spm.spatial.preproc.tissue(5).warped = [0 1];
matlabbatch{1}.spm.spatial.preproc.tissue(6).tpm = {'/oscar/data/ofeldman/alamba1/spm12/tpm/TPM.nii,6'};
matlabbatch{1}.spm.spatial.preproc.tissue(6).ngaus = 2;
matlabbatch{1}.spm.spatial.preproc.tissue(6).native = [0 0];
matlabbatch{1}.spm.spatial.preproc.tissue(6).warped = [0 1];
matlabbatch{1}.spm.spatial.preproc.warp.mrf = 1;
matlabbatch{1}.spm.spatial.preproc.warp.cleanup = 1;
matlabbatch{1}.spm.spatial.preproc.warp.reg = [0 0.001 0.5 0.05 0.2];
matlabbatch{1}.spm.spatial.preproc.warp.affreg = 'mni';
matlabbatch{1}.spm.spatial.preproc.warp.fwhm = 0;
matlabbatch{1}.spm.spatial.preproc.warp.samp = 3;
matlabbatch{1}.spm.spatial.preproc.warp.write = [0 0];

% Run the batch
spm_jobman('run', matlabbatch);
clear matlabbatch


