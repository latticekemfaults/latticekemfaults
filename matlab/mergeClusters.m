function cout = mergeClusters(c1, c1idx, c2, c2idx, l)
%MERGECLUSTERS Summary of this function goes here
%   Detailed explanation goes here

%perform marginalization, if needed
cout = c1;
if exist('c1idx', 'var') && not(isempty(c1idx))
    [cout.kg, cout.pkg] = marginalizeCluster(c1, c1idx);
    cout.idx = cout.idx(c1idx);
end

cm = c2;
if exist('c2idx', 'var') && not(isempty(c2idx))
    [cm.kg, cm.pkg] = marginalizeCluster(c2, c2idx);
    cm.idx = cm.idx(c2idx);
end

% merge indizes
cout.idx = [cout.idx, cm.idx];

% merge key guesses and probabilites
[cout.kg, cout.pkg] = generateKeyguesses({cout.kg, cm.kg}, {cout.pkg, cm.pkg});

% set new l
if exist('l', 'var')
    cout.l = l(cout.idx, :);
end


end

