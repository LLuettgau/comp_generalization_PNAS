function all_choices = get_choices(pch, subpch, obs_choice, correct_answer, simulation)


    %% fitting parameters to data
    if simulation == 0
        %compare subject's choices with softmax generated choice probabilities
        %(pch) in a new vector subpch --> how likely is the model's
        %(Softmax) prediction of a choice given the observed
        %choices of a participant

        for k = 1:size(subpch,1)
            if obs_choice(k) == 1 %if subject chose left, take left softmax generated choice probability
                subpch(k,1) = pch(k,1);
            elseif obs_choice(k) == 2 %if subject chose right, take right softmax generated choice probability
                subpch(k,1) = pch(k,2);
            elseif isnan(obs_choice(k)) %if subject did not respond - NaN
                subpch(k,1) = NaN;
            end
    
            if subpch(k,1) == 0
                subpch(k,1) = 1e-8; %replace 0 probability by extremely small value (problematic to have 0 in sum of logs)
            end
        end
    
        chprob = pch;
        sim_choices = [];
        sim_responses_corr = [];
        
    %% simulations
    else

        rand_choices = rand(size(subpch,1),1);
        sim_choices = NaN(size(rand_choices));
        sim_choices(rand_choices<pch(:,1)) = 1;
        sim_choices(rand_choices>pch(:,1)) = 2;

        clear chosen_option_correct
        chosen_option_correct(:,1) = sim_choices == correct_answer;

        sim_responses_corr = double(chosen_option_correct(:,1));
        sim_responses_corr(find(isnan(sim_choices))) = NaN;

        chprob = pch;
        subpch = [];
    end

    %prepare output
    all_choices.chprob = chprob;
    all_choices.sim_responses_corr = sim_responses_corr;
    all_choices.sim_choices = sim_choices;
    all_choices.subpch = subpch;


