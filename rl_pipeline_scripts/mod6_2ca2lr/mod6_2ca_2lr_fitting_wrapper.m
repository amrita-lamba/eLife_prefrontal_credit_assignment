%% Model fitting wrapper for RL with logistic choice function 
% model uses two additional parameters for credit assignment and the learning rate  

% read in the data 
data = readtable('adult_behavData.csv');

% initalize random seed 
rng shuffle; 
params = struct;

%lower bounds, vector with one entry per param 
% LRPos, LRNeg, bias, m, prior, caPos,caNeg 
lb = [0,0,0,0.2,0,0,0]; 
ub = [1,1,2,2,2,1,1]; 

% initialize structure to save best-fitting parameters
params.subjNum = [];
params.subjCond = []; 
params.subjLRPos = []; 
params.subjLRNeg = []; 
params.subjBias = [];
params.subjM = []; 
params.subjCaPos = [];
params.subjCaNeg = [];
params.subjPrior = []; 
params.subjBayesLL = []; 
params.subjAIC = [];

% save deltas and q-vals for each subject
deltas_allSubs = zeros(60,60); 
qvals_allSubs = zeros(60,60); 

sub_list = unique(data.subject);

for cond = 1:2

    for sub = 1:length(sub_list)

        subjData = data(data.subject == sub & data.condition == cond,:);       

        % set number of iterations
        nIter = 40; 
        rng shuffle; 

        % create array to save trial-wise deltas and q-values from model 
        deltas_allIterations = zeros(60,nIter); 
        qval_allIterations = zeros(60,nIter);

        for iter = 1:nIter   

            % randomly select starting points for model within parameter bounds 
            initVals = [rand rand rand*2 0.2+rand rand*2 rand rand];
            [res(iter,:),lik(iter),flag,out,lambda,grad,hess] = ...
                fmincon(@(x) mod6_2ca_2lr_loglikelihood(subjData,x(1),x(2),x(3),x(4),x(5),x(6),x(7)), initVals,[],[],[],[],lb,ub);

            % load in deltas and qvals from iteration
            load('deltas.mat'); 
            load('q_vals.mat'); 
            deltas_allIterations(:,iter) = delta_vec; 
            qval_allIterations(:,iter) = Q_vec;  

        end  

        minBayesLL = min(lik);
        extractParams = find(lik == minBayesLL);

        % if more than 1 equivalent sets with min BayesLL take first
        if length(extractParams) > 1
            extractParams = extractParams(1);
        end 

        % save best fitting deltas for each subject
        deltas_allSubs(:,sub) = deltas_allIterations(:,extractParams); 
        qvals_allSubs(:,sub)  = qval_allIterations(:,extractParams); 

        LRPos = res(extractParams,1); 
        LRNeg = res(extractParams,2); 
        bias  = res(extractParams,3);
        m     = res(extractParams,4); 
        prior = res(extractParams,5); 
        caPos = res(extractParams,6); 
        caNeg = res(extractParams,7); 

        numParams = length(lb);

        % calculate the AIC using BayesLL. NOTE: This script calcualtes AIC
        % in negative space so penalization term is subtracted. 
        % Use +2*numParams if LL is negative 
        AIC_val = -2*(minBayesLL)-2*numParams;

        % save parameters 
        params.subjNum     = [params.subjNum sub]; 
        params.subjCond    = [params.subjCond cond]; 
        params.subjLRPos   = [params.subjLRPos LRPos]; 
        params.subjLRNeg   = [params.subjLRNeg LRNeg]; 
        params.subjBias    = [params.subjBias bias];
        params.subjM       = [params.subjM m];
        params.subjPrior   = [params.subjPrior prior];
        params.subjCaPos   = [params.subjCaPos caPos];
        params.subjCaNeg   = [params.subjCaNeg caNeg];
        params.subjBayesLL = [params.subjBayesLL minBayesLL];
        params.subjAIC     = [params.subjAIC AIC_val];

        save(['mod6_2ca_2lr.mat'],'params');

    end

end 

% create data table 
t = [params.subjNum' params.subjCond' params.subjLRPos' params.subjLRNeg' params.subjBias' params.subjM' params.subjPrior' params.subjCaPos' params.subjCaNeg' params.subjBayesLL' params.subjAIC'];
t_table = array2table(t);

% add table headers 
t_table.Properties.VariableNames = {'subject','condition','LRPos','LRNeg','bias','m','prior','caPos','caNeg','BayesLL','AIC'};

filename = 'mod6_2ca_2lr.csv';                                   
writetable(t_table,filename,'WriteVariableNames',true);   

%% save deltas and model q-vals

num_trials = 60;
num_subs   = 30; 

trial_col = repmat([1:num_trials],1,num_trials)';
sub_col   = [(repelem([1:num_subs],1,num_trials))';(repelem([1:num_subs],1,num_trials))'];
cond_col  = [(repelem([1],1,num_trials*num_subs))';(repelem([2],1,num_trials*num_subs))'];
delta_col = []; % store model delta
Q_col     = []; % store q-val

for cond = 1:2
    for sub = 1:length(sub_list)

        curDeltas = deltas_allSubs(:,sub); 
        delta_col = [delta_col; curDeltas]; 

        curQ = qvals_allSubs(:,sub); 
        Q_col = [Q_col; curQ];

    end 
end 

delta_q_table = [sub_col, cond_col, trial_col, delta_col,Q_col]; 
delta_q_table = array2table(delta_q_table);  

% add table headers 
delta_q_table.Properties.VariableNames = {'subject','condition','trial','delta','q_val'};

filename = 'mod6_2ca_2lr_deltas_q.csv';                                   
writetable(delta_q_table,filename,'WriteVariableNames',true);   

