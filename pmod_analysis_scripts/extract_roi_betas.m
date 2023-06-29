function [individual_b,pc] = extract_roi_betas(R, design)

    % fetch data
    Y = get_marsy(R, design, 'mean');

    % get contrasts from original design 
    xCon = get_contrasts(design); 

    % estimate design on ROI
    E = estimate(design, Y); 

    % put contrasts from original design back into design object
    E = set_contrasts(E, xCon); 

    % get design betas 
    b = betas(E); 

    % acquire statistics 
    marsS = compute_contrasts(E, 1:length(xCon)); 

    % individual betas
    individual_b = marsS.MVres.y_obs; 
    
    % corrected pvalue (corrected for number of ROIs)
    pc = marsS.Pc(1); 
    
%     % get time course 
%     mY = get_marsy(roi, data, 'mean'); % extract data into marsy data object
%     y = summary_data(mY); % get summary time course(s)
    
end