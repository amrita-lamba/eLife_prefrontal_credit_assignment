function al_proproc_RSA(subNum,scriptPath)

% Default directories
subjPath = [scriptPath,'/NIFTI_2mm/',int2str(subNum)];

rundirs = dir(fullfile(subjPath, '*ep2d_functional*'));
rundirs = struct2cell(rundirs);
rundirs = rundirs(1,:);

mpragedir =  dir(fullfile(subjPath, '*mprage*'));
mpragedir = struct2cell(mpragedir);
mpragedir = mpragedir(1,:);

% **********************************************************************

cd(subjPath)

%% Slice Time Correction 

% find scans for each run
data = {};
for ri = 1:4

runPath = fullfile(rundirs{ri});
scanfnames = dir(fullfile(subjPath,runPath, '4*.img'));

data = [data; ...
    reshape(strcat(runPath, '/', {scanfnames.name}, ',1'), [], 1)];    
    
end

matlabbatch = {};
matlabbatch{1}.spm.temporal.st.scans = {data}';
matlabbatch{1}.spm.temporal.st.nslices = 38;
matlabbatch{1}.spm.temporal.st.tr = 2;
matlabbatch{1}.spm.temporal.st.ta = 1.94736842105263;
matlabbatch{1}.spm.temporal.st.so = [1 3 5 7 9 11 13 15 17 19 21 23 25 27 29 31 33 35 37 2 4 6 8 10 12 14 16 18 20 22 24 26 28 30 32 34 36 38];
matlabbatch{1}.spm.temporal.st.refslice = 1;
matlabbatch{1}.spm.temporal.st.prefix = 'a';


% Run the batch
spm_jobman('run', matlabbatch);
clear matlabbatch


%% Realignment

% find scans for each run
data = {};
for ri = 1:4

runPath = fullfile(rundirs{ri});
scanfnames = dir(fullfile(subjPath,runPath, 'a4*.img'));

data = [data; ...
    reshape(strcat(runPath, '/', {scanfnames.name}, ',1'), [], 1)];

end 

matlabbatch = {};
matlabbatch{1}.spm.spatial.realign.estwrite.data = {data}'; 
matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.quality = 0.9;
matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.sep = 3;
matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.fwhm = 5;
matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.rtm = 1;
matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.interp = 2;
matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.wrap = [0 0 0];
matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.weight = '';
matlabbatch{1}.spm.spatial.realign.estwrite.roptions.which = [2 1];
matlabbatch{1}.spm.spatial.realign.estwrite.roptions.interp = 4;
matlabbatch{1}.spm.spatial.realign.estwrite.roptions.wrap = [0 0 0];
matlabbatch{1}.spm.spatial.realign.estwrite.roptions.mask = 1;
matlabbatch{1}.spm.spatial.realign.estwrite.roptions.prefix = 'r';

% Run the batch
spm_jobman('run', matlabbatch);
clear matlabbatch


%% Coregistration

% find reference image 
runPath = fullfile(rundirs{1});
scanfnames = dir(fullfile(subjPath,runPath, 'meana*.img'));
refImage = {}; 
refImage = [refImage; ...
    reshape(strcat(runPath, '/', {scanfnames.name}, ',1'), [], 1)];    

% find source image 
anatPath = fullfile(mpragedir{1});
scanfnames = dir(fullfile(subjPath,anatPath, '4*.img'));
sourceImage = {}; 
sourceImage = [sourceImage; ...
    reshape(strcat(anatPath, '/', {scanfnames.name}), [], 1)];  

matlabbatch = {};
matlabbatch{1}.spm.spatial.coreg.estwrite.ref = refImage; 
matlabbatch{1}.spm.spatial.coreg.estwrite.source = sourceImage;
matlabbatch{1}.spm.spatial.coreg.estwrite.other = {''};
matlabbatch{1}.spm.spatial.coreg.estwrite.eoptions.cost_fun = 'nmi';
matlabbatch{1}.spm.spatial.coreg.estwrite.eoptions.sep = [4 2];
matlabbatch{1}.spm.spatial.coreg.estwrite.eoptions.tol = [0.02 0.02 0.02 0.001 0.001 0.001 0.01 0.01 0.01 0.001 0.001 0.001];
matlabbatch{1}.spm.spatial.coreg.estwrite.eoptions.fwhm = [7 7];
matlabbatch{1}.spm.spatial.coreg.estwrite.roptions.interp = 4;
matlabbatch{1}.spm.spatial.coreg.estwrite.roptions.wrap = [0 0 0];
matlabbatch{1}.spm.spatial.coreg.estwrite.roptions.mask = 0;
matlabbatch{1}.spm.spatial.coreg.estwrite.roptions.prefix = 'r';

% Run the batch
spm_jobman('run', matlabbatch);
clear matlabbatch

%% Segmentation 

% find segmentation volume
anatPath = fullfile(mpragedir{1});
scanfnames = dir(fullfile(subjPath,anatPath, 'r4*.img'));
volImage = {}; 
volImage = [volImage; ...
    reshape(strcat(anatPath, '/', {scanfnames.name}, ',1'), [], 1)];  

matlabbatch = {};
matlabbatch{1}.spm.spatial.preproc.channel.vols = volImage; 
matlabbatch{1}.spm.spatial.preproc.channel.biasreg = 0.001;
matlabbatch{1}.spm.spatial.preproc.channel.biasfwhm = 60;
matlabbatch{1}.spm.spatial.preproc.channel.write = [0 1];
matlabbatch{1}.spm.spatial.preproc.tissue(1).tpm = {'/Users/astrocyte/spm12/tpm/TPM.nii,1'};
matlabbatch{1}.spm.spatial.preproc.tissue(1).ngaus = 1;
matlabbatch{1}.spm.spatial.preproc.tissue(1).native = [1 0];
matlabbatch{1}.spm.spatial.preproc.tissue(1).warped = [0 0];
matlabbatch{1}.spm.spatial.preproc.tissue(2).tpm = {'/oscar/data/ofeldman/alamba1/spm12/tpm/TPM.nii,2'};
matlabbatch{1}.spm.spatial.preproc.tissue(2).ngaus = 1;
matlabbatch{1}.spm.spatial.preproc.tissue(2).native = [1 0];
matlabbatch{1}.spm.spatial.preproc.tissue(2).warped = [0 0];
matlabbatch{1}.spm.spatial.preproc.tissue(3).tpm = {'/oscar/data/ofeldman/alamba1/spm12/tpm/TPM.nii,3'};
matlabbatch{1}.spm.spatial.preproc.tissue(3).ngaus = 2;
matlabbatch{1}.spm.spatial.preproc.tissue(3).native = [1 0];
matlabbatch{1}.spm.spatial.preproc.tissue(3).warped = [0 0];
matlabbatch{1}.spm.spatial.preproc.tissue(4).tpm = {'/oscar/data/ofeldman/alamba1/spm12/tpm/TPM.nii,4'};
matlabbatch{1}.spm.spatial.preproc.tissue(4).ngaus = 3;
matlabbatch{1}.spm.spatial.preproc.tissue(4).native = [1 0];
matlabbatch{1}.spm.spatial.preproc.tissue(4).warped = [0 0];
matlabbatch{1}.spm.spatial.preproc.tissue(5).tpm = {'/oscar/data/ofeldman/alamba1/spm12/tpm/TPM.nii,5'};
matlabbatch{1}.spm.spatial.preproc.tissue(5).ngaus = 4;
matlabbatch{1}.spm.spatial.preproc.tissue(5).native = [1 0];
matlabbatch{1}.spm.spatial.preproc.tissue(5).warped = [0 0];
matlabbatch{1}.spm.spatial.preproc.tissue(6).tpm = {'/oscar/data/ofeldman/alamba1/spm12/tpm/TPM.nii,6'};
matlabbatch{1}.spm.spatial.preproc.tissue(6).ngaus = 2;
matlabbatch{1}.spm.spatial.preproc.tissue(6).native = [0 0];
matlabbatch{1}.spm.spatial.preproc.tissue(6).warped = [0 0];
matlabbatch{1}.spm.spatial.preproc.warp.mrf = 1;
matlabbatch{1}.spm.spatial.preproc.warp.cleanup = 1;
matlabbatch{1}.spm.spatial.preproc.warp.reg = [0 0.001 0.5 0.05 0.2];
matlabbatch{1}.spm.spatial.preproc.warp.affreg = 'mni';
matlabbatch{1}.spm.spatial.preproc.warp.fwhm = 0;
matlabbatch{1}.spm.spatial.preproc.warp.samp = 3;
matlabbatch{1}.spm.spatial.preproc.warp.write = [0 1];

% Run the batch
spm_jobman('run', matlabbatch);
clear matlabbatch

%% Normalization

% find deformation field image
anatPath = fullfile(mpragedir{1});
scanfnames = dir(fullfile(subjPath,anatPath, 'y_r4*.nii'));
deformationImage = {}; 
deformationImage = [deformationImage; ...
    reshape(strcat(anatPath, '/', {scanfnames.name}), [], 1)];  

% find scans to normalize
data = {};
for ri = 1:4
    
runPath = fullfile(rundirs{ri});
scanfnames = dir(fullfile(subjPath,runPath, 'ra4*.img'));

data = [data; ...
    reshape(strcat(runPath, '/', {scanfnames.name}, ',1'), [], 1)];    
    
end

matlabbatch = {};
matlabbatch{1}.spm.spatial.normalise.write.subj.def = deformationImage; 
matlabbatch{1}.spm.spatial.normalise.write.subj.resample = data;
matlabbatch{1}.spm.spatial.normalise.write.woptions.bb = [-78 -112 -70
                                                          78 76 85];
matlabbatch{1}.spm.spatial.normalise.write.woptions.vox = [2 2 2];
matlabbatch{1}.spm.spatial.normalise.write.woptions.interp = 4;
matlabbatch{1}.spm.spatial.normalise.write.woptions.prefix = 'w';

% Run the batch
spm_jobman('run', matlabbatch);
clear matlabbatch

%% Smoothing

% find scans to smooth
data = {};
for ri = 1:4
    
runPath = fullfile(rundirs{ri});
scanfnames = dir(fullfile(subjPath,runPath, 'wra4*.img'));

data = [data; ...
    reshape(strcat(runPath, '/', {scanfnames.name}, ',1'), [], 1)];    
    
end

matlabbatch = {};
matlabbatch{1}.spm.spatial.smooth.data = data; 
matlabbatch{1}.spm.spatial.smooth.fwhm = [2 2 2]; % smoothing kernal
matlabbatch{1}.spm.spatial.smooth.dtype = 0;
matlabbatch{1}.spm.spatial.smooth.im = 0;
matlabbatch{1}.spm.spatial.smooth.prefix = 's';

% Run the batch
spm_jobman('run', matlabbatch);
clear matlabbatch

