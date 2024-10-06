

function Aout = find_best_perm(A, B)

    maxsteps = 20; %20

    % define distance
    % use correlation
    distfn = @(A, B) corr(A(:), B(:));

    [n, ~] = size(B);
    out = cell(5, 1);
    for rep = 1:5 %5
        perm = randperm(n);
        orig_err = -Inf;
        best_err = -Inf;
        for step = 1:maxsteps
            best = [1, 1];
            orig_err = distfn(A(perm, perm), B);
            for i = 1:n
                for j = 1:i-1
                    perm_tmp = perm;
                    perm_tmp([i, j]) = perm_tmp([j, i]);
                    err = distfn(A(perm_tmp, perm_tmp), B);
                    
                    if (err > orig_err) && (err > best_err)
                    % if (err < orig_err) && (err < best_err)
                        best_err = err;
                        best = [i, j];
                    end
                end
            end
            i = best(1);
            j = best(2);
            if all(best == [1, 1])
                out{rep} = {A(perm, perm), best_err, perm, step};
                break;
            end
            perm([i, j]) = perm([j, i]);
        end
    end
    [~, i] = max(cellfun(@(o) o{2}, out));
    %[~, i] = min(cellfun(@(o) o{2}, out));

    Aout = out{i}{1};
end

% function Aout = find_best_perm(A, B)
% 
% %Kim's suggestions
% %learning phase
% % 1. learn SF (or SR) matrix M
% % SVD of M:
% % M = USV'
% %
% % transfer phase(s)
% % 2. learn new M" in transfer
% %
% % M" = US"V'
% % U, V from learning phase
% % S" = U' M" V
% %
% % off-diagonal values of S" constitute lack of transfer, aka variance of M" not captured by M's principal components
% % you can approximate transfer by keeping only diagonal entries of S" (variance that *is* captured by U, V)
% 
%     [U_prior,S_prior,V_prior] = svd(A);
% 
%     %calculate transfer phase singular value matrix, only
%     %retaining the diagonal values (i.e., what is shared!)
%     S_trans = U_prior' * B * V_prior;
%     S_trans_diag = diag(diag(S_trans));
% 
%     %new inferred transfer phase SF matrix
%     Aout = U_prior * S_trans_diag * V_prior';
% 
% 
% 
% end



