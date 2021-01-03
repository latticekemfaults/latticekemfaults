function h = entropy(probs)
%ENTROPY Shannon entropy of a probability distribution.
%   h = ENTROPY(probs) computes the Shannon entropy of the discrete 
%   probabilty distribution given by the probabilty vector probs.

probs = probs/sum(probs);
l = log2(probs);
l(probs == 0) = 0; % have to fix the behavior for zeros
h = -sum(probs.*l);

end

