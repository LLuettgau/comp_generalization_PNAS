function succ_prob = get_SF_probs(SF_M, probe_features, options, iTrial)

        succ_prob = [SF_M(probe_features(iTrial,1), options(iTrial,1)) + ...
                     SF_M(probe_features(iTrial,1), options(iTrial,2)) + ...
                     SF_M(probe_features(iTrial,2), options(iTrial,1)) + ...
                     SF_M(probe_features(iTrial,2), options(iTrial,2)) ...
                     SF_M(probe_features(iTrial,1), options(iTrial,3)) + ...
                     SF_M(probe_features(iTrial,1), options(iTrial,4)) + ...
                     SF_M(probe_features(iTrial,2), options(iTrial,3)) + ...
                     SF_M(probe_features(iTrial,2), options(iTrial,4))];



       % %left option
       % choice_opt_left = [SF_M(probe_features(iTrial,1),options(iTrial,[1 2])); ...
       %                    SF_M(probe_features(iTrial,2),options(iTrial,[1 2]))];
       % 
       % succ_prob(1) = sum(sum(choice_opt_left,2));
       % % % choice_opt_left = sum(choice_opt_left,1);
       % % % succ_prob(1) = norm(choice_opt_left);
       % % 
       % %right option
       % choice_opt_right = [SF_M(probe_features(iTrial,1),options(iTrial,[3 4])); ...
       %                     SF_M(probe_features(iTrial,2),options(iTrial,[3 4]))];
       % 
       % succ_prob(2) = sum(sum(choice_opt_right,2));
       % % choice_opt_right = sum(choice_opt_right,1);
       % % succ_prob(2) = norm(choice_opt_right);

end