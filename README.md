# eLife_fMRI_pipeline

This repository contains the behavioral and analyzed neural data for RSA and parametric modulation analyses reported in Lamba, Nassar, FeldmanHall, 2023, eLife.

Details about the data are provided below:

RSA_analyis_scripts (Fig. 4-6 of paper)
Behavioral data is stored in "adult_behavData_2ca2lr.csv"

subjNum          - initial ID assigned to subject during the scan session
condition        - 1 = Trust Task, 2 = Bandit Task
runNum           - specifies whether data was collected from first or second run during scan session
trial            - trial number in task (60 trials per task)
playerNum        - partner type: 1 = low return, 2 = high return, 3 = neutral, 4 = random
endowment        - subject monetary investment/choice on trial: NOTE: -1 values in endowment, returned, and trial RT columns denote that subject did not respond on that trial
propReturned     - preset return rate for partner type on trial. These were randomly shuffled and assigned to participants in task code
returned         - monetary amount returned to subject on trial during feedback
trialRT          - response time during choice
optimal_choice   - payoff maximizing response on current trial
binned_trial_num - trial number with specific partner type
delta            - prediction error on current trial from 2ca2lr model
q-val            - Q-value/V on current trial from 2ca2lr model

Details about subjects for analysis purposes stored in "Adult_Run_List.csv"
counterBalance      - 1 = Trust Task completed first,  2 = Bandit Task completed first
exclusionCode       - 1 = subject excluded from final analyses because head motion in one direction exceeded 3mm cutoff, 0 = subject included in final analyses
runxrun_correction  - 1 = realignment to first scan of every run if motion/body adjustment is excessive between runs: Note, this was not needed for any subjects in sample

Results from RSA choice/feedback phase analyses:
"choice_identity_rsa_roi.csv" & "feedback_identity_rsa_roi.csv"
model  - predictor RDM
roi    - region of interest (see ROI_coords.xlsx for size of clusters and location)
b_est  - beta coefficient from RSA
t_stat - t-statistic from RSA

Results from cross-timepoint analyses
"alt_conj_rho.xlsx"
roi - conjunction ROI
h - predictor RDM denoting whether RSA was computed from choice, feedback, or across timepoints using even and odd trials (crossphase)
rho - correlation coefficient from RSA

Timings:
All onsets, offsets, and durations for the task are stored in the Full_DM_Timings_Combined folder. Folder also contains beta images from deconvolution and the first-level design matrix (SPM.mat)

pmod_analysis_scripts (Fig. 7 of paper)
Contains the same behavioral and run list files as RSA folder
Results of parametric modulation analyses: "delta_beta_table.xlsx"
beta         - observed magnitude of BOLD response for subject in ROI, using marsbar toolbox for estimation
p_corrected  - adjusted p-value from marsbar toolbox, p-value corrected for number of ROIs
brain_region - ROI
contrast     - species which condition results are from, TG = Trust Task, SM = Bandit task
