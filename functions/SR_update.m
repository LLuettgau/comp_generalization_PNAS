function [M, SF_M, sim_res, inf_sim_res, sim_choices, inf_sim_choices] = SR_update(sched, alpha, gamma, tau, w, omega, dims_SF, choice_trials, mechanism, M_prior, condition)

    if mechanism < 3 || mechanism == 4
        M = zeros(dims_SF);
    else
        M = zeros(dims_SF.SR);
        SF_M = zeros(dims_SF.SF);
        sched_tmp = sched;
        sched = sched.SR;
    end
    
    iTrial = 1;
    k = zeros(size(M,1),1);
    state_visit_counts = zeros(dims_SF);
    for t=1:(size(sched,1)-1)
        M(:,:,t+1)= M(:,:,t);
        
        if mechanism == 3
            SF_M(:,:,t+1)= SF_M(:,:,t);
        end
 

        if mechanism == 1 %SF compound updating

            onehot = zeros(1,dims_SF(1));
            onehot(sched(t)) = 1;

            M(sched(t),:,t+1)= M(sched(t),:,t)+ alpha*[onehot+ gamma*M(sched(t+1),:,t) - M(sched(t),:,t)];


        elseif mechanism == 0 || mechanism == 2 || mechanism == 4 %SF feature updating - reusing prior learning phase SR or vanilla

            twohot = zeros(1,dims_SF(1));
            twohot(sched(t,:)) = 1;

            for ii=1:2
                M(sched(t,ii),:,t+1)= M(sched(t,ii),:,t) + .5*alpha*[twohot+ gamma*sum(M(sched(t+1,:),:,t),1) - sum(M(sched(t,:),:,t),1)];

                %count of state visitations, for exponential decay
                k(sched(t,ii)) = k(sched(t,ii)) + 1;
            end

        elseif mechanism == 3 %SF compound/SF feature hybrid updating

            sched_SR = sched_tmp.SR;
            sched_SF = sched_tmp.SF;
            
            %SF compound update
            onehot_SR = zeros(1,dims_SF.SR(1));
            onehot_SR(sched_SR(t)) = 1;

            M(sched_SR(t),:,t+1)= M(sched_SR(t),:,t)+ alpha*[onehot_SR+ gamma*M(sched_SR(t+1),:,t) - M(sched_SR(t),:,t)];

            %SF feature update
            twohot = zeros(1,dims_SF.SF(1));
            twohot(sched_SF(t,:)) = 1;

            for ii=1:2
                SF_M(sched_SF(t,ii),:,t+1)= SF_M(sched_SF(t,ii),:,t)+ .5*alpha*[twohot+ gamma*sum(SF_M(sched_SF(t+1,:),:,t),1) - sum(SF_M(sched_SF(t,:),:,t),1)];
            end
            
        end


        %make choices
        if mod(t,12) == 0 %probes occur after every 12th trial
            
            if mechanism < 3
                SF_M = M;
            end

            %update M reusing prior learning phase SR
            if mechanism == 0 && ~isempty(M_prior)
                M_trans = M(:,:,end);
                %find permutation that maps prior learning SR to new transfer
                %learning SR

%                 % toy example
%                 B = [1 0 0 0; 0 1 0 0; 0 0 1 0; 0 0 0 1];
                % B = [0 1 0 0; 0 0 1 0; 0 0 0 1; 1 0 0 0];
                % A = [0 0 0 1; 1 0 0 0; 0 1 0 0; 0 0 1 0];
                % B = [0 1 0 0 0 0; 0 0 1 0 0 0; 0 0 0 1 0 0; 0 0 0 0 1 0; 0 0 0 0 0 1; 1 0 0 0 0 0];
                % A = [0 0 0 0 0 1; 1 0 0 0 0 0; 0 1 0 0 0 0; 0 0 1 0 0 0; 0 0 0 1 0 0; 0 0 0 0 1 0];
                % 
                % B = B * 5;
                % A = A * 5;
                % 
                % B_SR = expm(B);
                % A_SR = expm(A);
                % 
                % Aout = find_optimal_perm(A, B);
                % A_SR_out = find_optimal_perm(A_SR, B_SR);

                %create transition matrix from SR                
                %T = (M^-1 - I) / -gamma
                %T = (M^-1 - I) / -.9;
                %get_TM = @(M) ((inv(M) - eye(size(M))) / -.9);

                % M_prior_permuted = find_best_perm(M_prior, M_trans);
                M_prior_permuted = find_optimal_perm(M_prior, M_trans, 20);

                % [A,B,R,U,V] = canoncorr(M_prior,M_trans);

                % %update transfer learning SR matrix
                % %same learning rate and mixing parameter
                % M_new = [];
                % for s = 1:size(M_trans,2)
                %     M_new(s,:) = (1 - alpha)^k(s) * M_prior_permuted(s,:) + (1-(1-alpha)^k(s)) * M_trans(s,:);
                % end
                
                % %using independent learning rate and mixing parameter mixing parameter
                % M_new = [];
                % for s = 1:size(M_trans,2)
                %     M_new(s,:) = exp(-w * k(s)) * M_prior_permuted(s,:) + (1 - exp(-w * k(s))) * M_trans(s,:);
                % end
                % 
                % %using fixed mixing parameter
                % for s = 1:size(M_trans,2)
                %     M_new(s,:) = (1 - 0.5)^k(s) * M_prior_permuted(s,:) + (1-(1-0.5)^k(s)) * M_trans(s,:);
                % end

                %using independent learning rate and static mixing parameter
                M_new = [];
                for s = 1:size(M_trans,2)
                    M_new(s,:) = w * M_prior_permuted(s,:) + (1 - w) * M_trans(s,:);
                end

                %weighted matrix becomes new transfer learning SR matrix 
                M(:,:,end) = M_new;

                % figure
                % subplot(1,4,1)
                % imagesc(M_prior); colorbar
                % 
                % subplot(1,4,2)
                % imagesc(M_prior_permuted); colorbar
                % 
                % subplot(1,4,3)
                % imagesc(M_trans); colorbar
                % 
                % subplot(1,4,4)
                % imagesc(M_new); colorbar
            elseif mechanism == 4 && ~isempty(M_prior)
                % %Tryout mixing in actual TM of transfer learning phase to validate
                % if iTrial == 1
                %     if condition == 1
                % 
                %         for i = 1:size(sched,1)-1
                % 
                %             current_state = sched(i,1);
                %             next_state = sched(i+1,1);
                % 
                %             %keep count of specific transition
                %             state_visit_counts(current_state,next_state) = state_visit_counts(current_state,next_state) + 1;
                %         end
                % 
                %     else
                % 
                %         for i = 1:size(sched,1)-1
                % 
                %             current_state = sched(i,2);
                %             next_state = sched(i+1,2);
                % 
                %             %keep count of specific transition
                %             state_visit_counts(current_state,next_state) = state_visit_counts(current_state,next_state) + 1;
                %         end
                % 
                %     end
                % 
                %     %make TM
                %     TM = state_visit_counts ./ max(state_visit_counts,1);
                % end
                % 
                % SF_M = TM;

                %usual case
                M_trans = M(:,:,end);
                SF_M = find_optimal_perm(M_prior, M_trans, 20);

            end
         
            %experience probes
            [sim_res_tmp, sim_choices_tmp, w] = make_choices(M(:,:,end), SF_M(:,:,end), tau, w, omega, choice_trials, mechanism, 1, iTrial);
            sim_res(iTrial,1) = sim_res_tmp;
            sim_choices(iTrial,1) = sim_choices_tmp;

            %inference probes
            [inf_sim_res_tmp, inf_sim_choices_tmp, w]  = make_choices(M(:,:,end), SF_M(:,:,end), tau, w, omega, choice_trials, mechanism, 2, iTrial);
            inf_sim_res(iTrial,1) = inf_sim_res_tmp;
            inf_sim_choices(iTrial,1) = inf_sim_choices_tmp;

            iTrial = iTrial + 1;

        end
        
    end

    if mechanism < 3
        SF_M = M;
    end
    
    %last trial choices
    %experience probes
    [sim_res_tmp, sim_choices_tmp, ~] = make_choices(M(:,:,end), SF_M(:,:,end), tau, w, omega, choice_trials, mechanism, 1, iTrial);
    sim_res(iTrial,1) = sim_res_tmp;
    sim_choices(iTrial,1) = sim_choices_tmp;

    %inference probes
    [inf_sim_res_tmp, inf_sim_choices_tmp, ~] = make_choices(M(:,:,end), SF_M(:,:,end), tau, w, omega, choice_trials, mechanism, 2, iTrial);
    inf_sim_res(iTrial,1) = inf_sim_res_tmp;
    inf_sim_choices(iTrial,1) = inf_sim_choices_tmp;


end