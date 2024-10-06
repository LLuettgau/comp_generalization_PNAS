function [Aout, optimalPerm, maxDist] = find_optimal_perm_old(A, B, nstarts)
    if nargin < 3
        nstarts = 20;
    end

    distfn = @(A, B) corr(A(:), B(:));

    perms = arrayfun(@(rep) find_perm(A, B), 1:nstarts, 'UniformOutput', false);
    dists = arrayfun(@(i) distfn(A(perms{i}, perms{i}), B), 1:nstarts);
    [maxDist, idx] = max(dists);
    optimalPerm = perms{idx};

    Aout = A(optimalPerm,optimalPerm);

end


function result = try_swap(A, B, perm, i, j, distfn)
    % create candidate permutation
    candperm = perm;
    candperm([i, j]) = candperm([j, i]);
    % return (inverse) distance
    result = distfn(A(candperm, candperm), B);
end

function perm = find_perm(A, B, maxsteps)
    if nargin < 3
        maxsteps = 5;
    end

    distfn = @(A, B) corr(A(:), B(:));
    [n, m] = size(B);
    perm = randperm(n);
    swaps = zeros(n*(n-1)/2, 2);
    count = 1;
    for i = 1:n
        for j = 1:i-1
            swaps(count, :) = [i, j];
            count = count + 1;
        end
    end

    for step = 1:maxsteps
        out = arrayfun(@(k) try_swap(A, B, perm, swaps(k, 1), swaps(k, 2), distfn), 1:size(swaps, 1));
        [~, idx] = max(out);
        j = swaps(idx,:);
        perm([j(1), j(2)]) = perm([j(2), j(1)]);

        % [i, j] = ind2sub(size(swaps), idx);
        % perm([i, j]) = perm([j, i]);
    end
end
