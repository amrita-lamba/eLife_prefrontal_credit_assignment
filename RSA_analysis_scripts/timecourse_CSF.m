function timecourse_CSF(subjnum,scriptPath)

% Default directories
subjPath = fullfile(scriptPath,'NIFTI/', int2str(subjnum));

% set up functional directories
rundirs = dir(fullfile(subjPath, '*ep2d_functional*'));
rundirs = struct2cell(rundirs);
rundirs = rundirs(1,:);

% set up anatomical directories 
mpragedir =  dir(fullfile(subjPath, '*mprage*'));
mpragedir = struct2cell(mpragedir);
mpragedir = mpragedir(1,:);

% **********************************************************************

cd(subjPath)
anatPath = fullfile(mpragedir{1});

%% Create binary ROI image from CSF map

% to create images use marsbar GUI to build ROIs for each subject (has to be done manually)
% load ROI image

cd(anatPath);
roi_file = dir(fullfile(subjPath,anatPath, '*roi.mat'));
load(roi_file.name);

%% extract all smoothed functional images

% navigate back to subj directory 
cd(subjPath);

% find scans to smooth
data = {};
for ri = 1:4
    
runPath = fullfile(rundirs{ri});
scanfnames = dir(fullfile(subjPath,runPath, 'swar4*.img'));

data = [data; ...
    reshape(strcat(runPath, '/', {scanfnames.name}), [], 1)];    
    
end

data = char(data); 

%% extract time course with marsbar 

mY = get_marsy(roi, data, 'mean'); % extract data into marsy data object
y = summary_data(mY); % get summary time course(s)

% save time course parameters
dlmwrite('CSF_timecourse.txt',y);


