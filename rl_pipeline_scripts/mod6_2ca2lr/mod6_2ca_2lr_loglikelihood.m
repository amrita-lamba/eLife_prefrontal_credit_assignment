%% loglikelihood function for 2ca2lr model  

function BayesLL = mod6_2ca_2lr_loglikelihood(subjData,LRPos,LRNeg,bias,m,prior,caPos,caNeg)

% track number of times participant has seen each partner to update Q matrix 
p1_cnt = 0;
p2_cnt = 0;
p3_cnt = 0;
p4_cnt = 0;

% use to keep track of partners not engaged with on current trial
players = (1:4);

% initialize the Q-matrix
Q_mat = zeros(15,4);

% estimate prior to initialize Q-matrix
Q_mat(1,:) = prior; 

% set up array to save deltas and Q-vals
delta_vec = zeros(60,1); 
Q_vec = zeros(60,1); 

% iterate through trials 
for i = 1:60

choice = subjData.endowment(i);    

    if choice == 0
    choice_idx = 1;
    elseif choice == 2.5
    choice_idx = 2;
    elseif choice == 5.0
    choice_idx = 3;
    elseif choice == 7.5
    choice_idx = 4;
    elseif choice == 10
    choice_idx = 5;
    elseif choice == -1
    choice_idx = 0;
    end 

    % find current partner 
    curPlayer = subjData.playerNum(i);

        if curPlayer == 1
            p1_cnt = p1_cnt + 1; 
            trial = p1_cnt;
        elseif curPlayer == 2
            p2_cnt = p2_cnt + 1; 
            trial = p2_cnt;
        elseif curPlayer == 3
            p3_cnt = p3_cnt + 1; 
            trial = p3_cnt;
        else 
            p4_cnt = p4_cnt + 1; 
            trial = p4_cnt;
        end 

    % update counts array     
    counts = [p1_cnt, p2_cnt, p3_cnt, p4_cnt]; 
        
        if choice ~= -1 
        
            if trial == 1

            % make a choice on trial one 
            predictedInvestment = (10/(1 + exp(-m*(prior-bias))));
            V = Q_mat((trial),curPlayer); 
            
            else 

            predictedInvestment = (10/(1 + exp(-m*(Q_mat((trial-1),curPlayer)-bias))));
            V = Q_mat((trial-1),curPlayer); 

            end 

        % discretize possible investments:
        possibleInvestments=0:2.5:10;  

        % last argument fixes sigma of normal distribution to 1
        choiceP = normpdf(possibleInvestments, predictedInvestment, 1);
        choiceP = choiceP./sum(choiceP);

        % compute trial loglikelihood
        trialLogLike(i) = log(choiceP(choice_idx));

        % calculate reward (returned - invested)
        reward = (choice*4*subjData.propReturned(i)) - choice;

        % calculate delta
        delta = reward - V;
        delta_vec(i) = delta; 

        % update Qs
        if delta >= 0 % positive outcome trials 
            
        Q_mat(trial,curPlayer) = V + LRPos*delta*caPos; 
        Q_vec(i) = Q_mat(trial,curPlayer); 
        
        % find players not engaged with on current trial 
        curPlayerIdx = find(players ~= curPlayer);

            if trial > 1 
                % implement credit assignment to players not seen on
                % current trial
                if counts(curPlayerIdx(1)) > 0
                Q_mat(counts(curPlayerIdx(1)),curPlayerIdx(1)) = Q_mat(counts(curPlayerIdx(1)),curPlayerIdx(1)) + (LRPos*delta*(1-caPos))/3;  
                end 

                if counts(curPlayerIdx(2)) > 0
                Q_mat(counts(curPlayerIdx(2)),curPlayerIdx(2)) = Q_mat(counts(curPlayerIdx(2)),curPlayerIdx(2)) + (LRPos*delta*(1-caPos))/3;
                end 

                if  counts(curPlayerIdx(3)) > 0
                Q_mat(counts(curPlayerIdx(3)),curPlayerIdx(3)) = Q_mat(counts(curPlayerIdx(3)),curPlayerIdx(3)) + (LRPos*delta*(1-caPos))/3;
                end 

            end    
          
        else % negative outcome trials 
          
        Q_mat(trial,curPlayer) = V + LRNeg*delta*caNeg; 
        Q_vec(i) = Q_mat(trial,curPlayer);
        
        % find players not seen 
        curPlayerIdx = find(players ~= curPlayer);
        
           if trial > 1 
                % add credit to Qs values for players not seen
                if counts(curPlayerIdx(1)) > 0
                Q_mat(counts(curPlayerIdx(1)),curPlayerIdx(1)) = Q_mat(counts(curPlayerIdx(1)),curPlayerIdx(1)) + (LRNeg*delta*(1-caNeg))/3;  
                end 

                if counts(curPlayerIdx(2)) > 0
                Q_mat(counts(curPlayerIdx(2)),curPlayerIdx(2)) = Q_mat(counts(curPlayerIdx(2)),curPlayerIdx(2)) + (LRNeg*delta*(1-caNeg))/3;
                end 

                if  counts(curPlayerIdx(3)) > 0
                Q_mat(counts(curPlayerIdx(3)),curPlayerIdx(3)) = Q_mat(counts(curPlayerIdx(3)),curPlayerIdx(3)) + (LRNeg*delta*(1-caNeg))/3;
                end 

            end    
        
        end 

        else % if subject did not respond on a trial other than first trial (Q-val doesn't change if sub misses first trial)
    
            if trial > 1
            Q_mat(trial,curPlayer) = Q_mat((trial-1),curPlayer);    
            Q_vec(i) = Q_mat((trial-1),curPlayer); 
            end 
            
        end
    
end 

BayesLL = -sum(trialLogLike);

% save deltas and qs
filename1 = 'deltas.mat';
filename2 = 'q_vals.mat';
save(filename1,'delta_vec');
save(filename2,'Q_vec');



