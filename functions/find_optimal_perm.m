function perm = find_optimal_perm(A, B, maxsteps)
    if nargin < 3
        maxsteps = 1000;
    end

    distfn = @(A, B) sum(abs(A(:) - B(:)).^2);

    [n, ~] = size(B);
    perm = randperm(n);
    swaps = zeros(n * (n - 1) / 2, 2);
    k = 1;
    for i = 1:n
        for j = 1:i-1
            swaps(k, :) = [i, j];
            k = k + 1;
        end
    end
    
    ll = try_swap(A, B, perm, 1, 1);
    
    for step = 1:maxsteps
        T = exp((1 - (5 * step + 1) / maxsteps));
        
        % try random swap
        idx = randi(size(swaps, 1));
        swap = swaps(idx, :);
        i = swap(1);
        j = swap(2);
        
        ll_cand = try_swap(A, B, perm, i, j);
        
        % accept with some randomness
        p = exp(-(ll_cand - ll) / T);
        
        if rand() < p
            % apply swap
            ll = ll_cand;
            perm([i, j]) = perm([j, i]);
        end
    end
end


function result = try_swap(A, B, perm, i, j)
    
    distfn = @(A, B) sum(abs(A(:) - B(:)).^2);

    % create candidate permutation
    candperm = perm;
    candperm([i, j]) = candperm([j, i]);
    % return (inverse) distance
    result = distfn(A(candperm, candperm), B);
end
