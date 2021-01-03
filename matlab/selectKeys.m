function [used,perr] = selectKeys(probin,err_threshold)
%SELECTKEYS Summary of this function goes here
%   Detailed explanation goes here

    numpkg = sum(cellfun(@(x) numel(x), probin));
    perr = zeros(1, numpkg - numel(probin)); % perr: store error prob when excluding n coeffs
    used = cellfun(@(x) {true(size(x))}, probin);
    
    usedx = used;
    % keep track of the current minimum value for each cluster
    pmin = cellfun(@(x) min(x), probin);
    
    ps = ones(size(probin)); % current correct prob for each cluster, starts with all ones
    for k = 1:numel(perr)
        % we want to minimize the increase in error prob at each step
        %   exclude 1 tuple in 1 cluster, correct prob gets multiplied with pnew/pold --> maximize this ratio at each step
        %   exclude entry with lowest probability: pold = pnew + pmin --> pnew/pold = (pold - pmin)/pold = 1 - pmin/pold
        %   --> minimize pmin/pold
        % 
        %   most likely entry never gets eliminated, because pmin/pold = 1 = maximum value
        [~, cidx] = min(pmin(:)./ps(:));
        [~, ridx] = min(probin{cidx}(used{cidx}));
        v = 1:numel(used{cidx});
        v = v(used{cidx});
        ridx = v(ridx);
        
        ps(cidx) = ps(cidx) - probin{cidx}(ridx);
        
        probin{cidx}(ridx) = 0;
        used{cidx}(ridx) = false;
        pmin(cidx) = min(probin{cidx}(used{cidx}));
        
        perr(k) = 1 - prod(ps);
        
        if perr(k) < err_threshold
            usedx{cidx}(ridx) = false;
        end
    end
    used = usedx;
end

