function [kgs, probs] = generateKeyguesses(components, componentprobs)
%GENKG Summary of this function goes here
%   Detailed explanation goes here

numcoeffs = cellfun(@(x) size(x, 2), components);
kgeach = cellfun(@(x) size(x, 1), components);
totkgs = prod(kgeach);

kgs = zeros(totkgs, sum(numcoeffs));

idx = cell(size(components));
inp = cellfun(@(x) {1:size(x, 1)}, components);
[idx{:}] = ndgrid(inp{:});

curridx = 1;
for k = 1:numel(components)
    matidx = (1:numcoeffs(k)) + (curridx - 1);
    kgs(:, matidx) = components{k}(idx{k}, :);
    curridx = curridx + numcoeffs(k);
end

if ~exist('componentprobs', 'var')
    probs = [];
else
    probmat = zeros(totkgs, numel(components));
    for k = 1:numel(components)
        probmat(:, k) = componentprobs{k}(idx{k}(:));
    end
    probs = prod(probmat, 2);
end

end

