function [sample] = inversionsample(vals, cdf, m, n)
%INVERSIONSAMPLE Sample values from an arbitrary discrete distribution.
%   sample = INVERSIONSAMPLE(vals, cdf, m, n) samples values from a discrete
%   distribution with support over vals and a cumulative distribution function cdf.
%   Samples [m x n] indepent values (default 1x1).

if ~exist('m', 'var')
    m = 1;
end

if ~exist('n', 'var')
    n = 1;
end

sample = zeros(m, n);

r = rand(m, n);

for i = 1:m
    for j = 1:n
        %find first entry where r < pmf
        ix = find(r(i, j) < cdf, 1, 'first');
        %and return the value
        sample(i, j) = vals(ix);
    end
end

end

