function [kgs, pkg] = marginalizeCluster(cluster, idx)
%MARGINALIZECLUSTER Marginalize probabilities of cluster
%   [kgs, pkg] = MARGINALIZECLUSTER(cluster, idx)
%   cluster is a struct describing the cluster, idx describes which coefficients to keep
%   idx is either a logical or a numeric array

numidx = size(cluster.kg, 2);
if ~islogical(idx)
    tmp = idx;
    idx = false(1, numidx);
    idx(tmp) = true;
end

[kgs, ~, ic] = unique(cluster.kg(:, idx), 'stable', 'rows');
pkg = accumarray(ic, cluster.pkg);


end

