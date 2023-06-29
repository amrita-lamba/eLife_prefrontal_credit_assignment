%% this script will delete all the intermediate preprocessing files to save room on disk drive 
%(will keep orignal NIFTY files)
function delete_interm_files(subNum,scriptPath)

% cd to the Timings folder
niftyPath = fullfile(scriptPath,'/NIFTI_2mm/', int2str(subNum));

% functional 
functionalPath = dir(fullfile(niftyPath, '*ep2d_functional*'));
functionalPath = struct2cell(functionalPath);
functionalPath = functionalPath(1,:);

%% this will loop through each subject's functional runs and delete the 
% intermediate files
cd(niftyPath);

    for ri = 1:length(functionalPath)
        curPath = fullfile(functionalPath{ri});
        cd(curPath);
        delete('a*.img');
        delete('a*.hdr');
        delete('ra*.img');
        delete('ra*.hdr');
        delete('wra*.img');
        delete('wra*.hdr');
        cd(niftyPath);
    end
    
