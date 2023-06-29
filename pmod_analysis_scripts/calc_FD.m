%% Calculate and plot Framewise Displacment
function calc_FD(subNum,runxrun_correction, scriptPath)

% set directories 
niftyPath = [scriptPath,'/NIFTI',int2str(subNum)];

if runxrun_correction == 0
    
    realignPath = dir(fullfile(niftyPath, '*ep2d_functional*'));
    realignPath = struct2cell(realignPath);
    realignPath = realignPath(1,1);
    realignPath = cell2mat(realignPath);

    realignFile = dir(fullfile(niftyPath, realignPath, '*rp*.txt'));
    realignFile = struct2cell(realignFile);
    realignFile = realignFile(1,1);
    realignFile = cell2mat(realignFile);

    % need to navigate to the folder to load the realign file 
    cd(fullfile(niftyPath,realignPath))
    realignParams = load(realignFile);

else 

    rundirs = dir(fullfile(subjPath, '*ep2d_functional*'));
    rundirs = struct2cell(rundirs);
    rundirs = rundirs(1,:);

    cd(subjPath)

    realignParams = [];

    for ri = 1:4
        cd(rundirs{ri})
        
        realignPath = dir(fullfile(niftyPath, '*ep2d_functional*'));
        realignPath = struct2cell(realignPath);
        realignPath = realignPath(1,1);
        realignPath = cell2mat(realignPath);

        %find the motion param file 
        realignFile = dir(fullfile(subjPath, rundirs{ri}, '*rp*.txt'));
        realignFile = struct2cell(realignFile);
        realignFile = realignFile(1,1);
        realignFile = cell2mat(realignFile);
        runParams = load(realignFile);

        %concatonate motion params across runs into single arrays 
        realignParams = [realignParams; runParams]; 
        cd(subjPath)
    end 

    % need to navigate to the folder to load the realign file 
    cd(fullfile(niftyPath,realignPath))
    
end 

%% set FD parameters 

% Note: Following section adapted from Cyril Pernet - University of Edinburgh code on
% github: https://github.com/CPernet/spmup/blob/master/QA/spmup_FD.m 
radius  = 50;

% de-mean and detrend
motion = spm_detrend(realignParams,1); 
% angle in radian by average head size = displacement in mm
motion(:,[4 5 6]) = motion(:,[4 5 6]).*radius; 

%% compute FD
cd(scriptPath)
D            = diff(motion,1,1); % 1st order derivative
D            = [zeros(1,6); D];  % set first row to 0
FD           = sum(abs(D),2); % framewise displacement a la Powers
RMS          = sqrt(mean(detrend(D).^2,2)); % root mean square for each column a la Van Dijk
FD_outliers  = spmup_comp_robust_outliers(FD);
RMS_outliers = spmup_comp_robust_outliers(RMS);

% save fd
cd(fullfile(niftyPath,realignPath))
dlmwrite('FD_params.txt',FD);
% dlmwrite('FD_outliers.txt',FD_outliers);

%% calculate threshold 

threshold = 1.2; % set threshold here 
FD_threshold = zeros(length(FD),1);
FD_exclude = find(FD > threshold); 
FD_prevFrame = FD_exclude-1;
FD_tplus1 = FD_exclude + 1; 
FD_tplus2 = FD_exclude + 2; 

FD_threshold(FD_exclude) = 1;
FD_threshold(FD_prevFrame) = 1;
FD_threshold(FD_tplus1) = 1;
FD_threshold(FD_tplus2) = 1;

FD_threshold = logical(FD_threshold);
dlmwrite('FD_threshold.txt',FD_threshold);

%% create a nuissance regressor for the first scan and last 4 scans of every run

% run_list = zeros(246,1); 
% 
% % add nuissance regressor for vol 1
% run_list(1) = 1; 
% 
% % add in nuissance regressor last 4 volumes
% run_list(243:246) = 1;
% 
% endPoint_regressors = repmat(run_list,4,1); 
% dlmwrite('endPoint_regressors.txt',endPoint_regressors);

%% save FD plot 

% current = pwd;
% fig     = 'save';
% 
% subplot(2,2,1); plot(motion(:,[1 2 3]),'LineWidth',3); axis tight; box on; grid on; title('translation')
% xlabel('Volumes'); ylabel('displacement in mm')
% subplot(2,2,2); plot(motion(:,[4 5 6]),'LineWidth',3); axis tight; box on; grid on; title('rotation')
% xlabel('Volumes'); ylabel('displacement in mm')
% subplot(2,2,3); plot(FD,'LineWidth',2); axis tight; box on; grid on; title('framewise displacement')
% hold on; tmp = FD_outliers.*FD; tmp(tmp==0)=NaN; plot(tmp,'ro','LineWidth',3);
% xlabel('Volumes'); ylabel('absolute displacement in mm')
% subplot(2,2,4); plot(RMS,'LineWidth',2); axis tight; box on; grid on; title('root mean squares')
% hold on; tmp = RMS_outliers.*RMS; tmp(tmp==0)=NaN; plot(tmp,'ro','LineWidth',3);
% xlabel('Volumes'); ylabel('average displacement in mm')
% 
% if strcmpi(fig,'save')
%     if exist(fullfile(niftyPath,'spm.ps'),'file')
%         print (gcf,'-dpsc2', '-bestfit', '-append', fullfile(niftyPath,'spm.ps'));
%     else
%         print (gcf,'-dpsc2', '-bestfit', '-append', fullfile(niftyPath,'spmup_QC.ps'));
%     end
%     close(gcf)
%     cd(current)
% end


