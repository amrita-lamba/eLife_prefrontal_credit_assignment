%% this script will delete all the interm files from preprocessing 
%(will keep orignal NIFTY files)
function delete_interm_files(subNum,scriptPath)

% cd to the Timings folder
niftyPath = fullfile(scriptPath,'/NIFTI_2mm/', int2str(subNum));

% functional 
functionalPath = dir(fullfile(niftyPath, '*ep2d_functional*'));
functionalPath = struct2cell(functionalPath);
functionalPath = functionalPath(1,:);

%% loop through each subject's functional run and delete the 
% intermediate files
cd(niftyPath);

    for ri = 1:length(functionalPath)
        curPath = fullfile(functionalPath{ri});
        cd(curPath);
        delete('r*.img');
        delete('r*.hdr');
        delete('ar*.img');
        delete('ar*.hdr');
        delete('war*.img');
        delete('war*.hdr');
        cd(niftyPath);
    end
    
