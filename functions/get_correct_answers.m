%compute correct responses for different probe types

function [correct_answer, probed_factor] = get_correct_answers(correct_answer, probed_factor, ...
                                              experimental_data, ...
                                                trial_type, ...
                                                probe_type, ...
                                                col_idx)
                                                                  
    %find probe trial idx
    trial_idx = find(strcmp(experimental_data.ScreenName, trial_type) ...
                             & strcmp(experimental_data.display, probe_type));


                          
    correct_answer(1:length(experimental_data.correctAnswer(trial_idx)),col_idx) = experimental_data.correctAnswer(trial_idx);



    %find probed factor
    probed_factor(1:length(experimental_data.correctAnswer(trial_idx)),col_idx) = experimental_data.featureProbe(trial_idx);
                                                                                                                                              
end



