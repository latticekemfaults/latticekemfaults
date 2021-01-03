%% choose parameter set and attack settings

 %#ok<*UNRCH>
 
kyber = 512; %this script was adapated from the attack on Kyber. For this reason, there are still some references to Kyber here
n = 2;

faultrange = 9000:500:19000;%newhope512
faultrange =  15000; %for demonstration
numfaults = 15000;
maxiterations = 10;

repeats = 1;

%exclusion and filtering stuff
err_threshold = 0.01;
domerge = true;
mergethreshold = 1000;

% errate = 0.01;
errate = 0;

plots = true;
store = false;

%%

addpath('ParforProgMon');

successes = false(numel(faultrange), repeats);
iterations = zeros(size(successes));
firsterror_best = zeros(size(successes));
firsterror_iterations = zeros(numel(faultrange), repeats, maxiterations);
sweep_repeat = 1;
sweep_numfaults = 1;

%% compute constants
% consts %initializes err_round_u and err_round_v

eta = 8;
vals = -eta:eta;

valsx = vals - vals(1);
pbino = binopdf(valsx, 2*eta, 0.5);

binoentropy = -sum(pbino.*log2(pbino));

%% sweep start
for sweep_repeat = 1:repeats
    t_repeat = tic;
    for sweep_numfaults = 1:numel(faultrange)
        numfaults = faultrange(sweep_numfaults);
        fprintf('SWEEP. NUMFAULTS: %d, REPEAT: %d\n', numfaults, sweep_repeat)
    

%% generate or read in test data
nx = numel(vals);

% kyber, simulated using reference implementation
f1 = 'ineqs.txt';
f2 = 'ineqs_es.txt';

% simulate
sstring = sprintf('./newhope%d_fault s -n %d', kyber, numfaults);
system(sstring);

movefile(f1, ['test_data/', f1]);
movefile(f2, ['test_data/', f2]);

% masked decoder, simulated
%     f1 = 'ineqs_masked.txt';
%     f2 = 'ineqs_es_masked.txt';

% ineqs using the fixed key
%     f1 = 'ineqs_fixed.txt';
%     f2 = 'fixed_es.txt';
%     end

d = readtable(['test_data/', f1], 'FileType', 'text', 'ReadVariableNames', false, 'MultipleDelimsAsOne', true, 'HeaderLines', 0);
l = single(table2array(d(:, 1:end-2)));
left = strcmp(table2cell(d(:, end-1)), '<')';
cutoff = single(table2array(d(:, end))');
cutoff = cutoff - 1; % correction mentioned above

numvars = size(l, 2);
fid = fopen(['test_data/', f2]);
es = textscan(fid, '%f', numvars);
es = single(es{1}');
fclose(fid);

numfaults = size(l, 1);

%introduce an error
if errate > 0
    x = binornd(1, errate, size(left)) > 0;
    left(x) = ~left(x);
end

%flip es --> se
nv2 = numvars/2;
se = [es(nv2+1:end), es(1:nv2)];
l = [l(:, nv2+1:end), l(:, 1:nv2)]';

% sanity check
realnoise = se*l;
checkleft = (realnoise <= cutoff);
check = isequal(checkleft, left);

if ~check
    fprintf('Data is inconsistent!\n'); % this is to be expected when using data from the real device
    err_rate = mean(checkleft ~= left);
    fprintf('Error rate: %.02f%%\n', err_rate*100);
end

 clear f1 f2 d numvars fid es realnoise checkleft check sstring


%% init data structure for key guesses
kgxcoeff = cell(2*kyber, 1);
for k = 1:numel(kgxcoeff)
    kgxcoeff{k} = vals';
end

kgx = cell(2*kyber/n, 1);

tmp = struct;
tmp.idx = 1:n;
tmp.kg = generateKeyguesses(kgxcoeff(1:n));
kgix = tmp.kg - vals(1) + 1;
tmp.pkg = prod(pbino(kgix), 2);
tmp.l = l(tmp.idx, :);

% r = n*eta*valssum(end);
r = n*eta*vals(end);
noiserange = -r:r;
tmp.dist = zeros(numel(noiserange), numfaults, 'single');

for k = 1:numel(kgx)
    kgx{k} = tmp;
    kgx{k}.idx = ((k-1)*n + 1):(k*n);
    kgx{k}.l = l(kgx{k}.idx, :);
%     kgx{k}.predshift = kgx{k}.kg*kgx{k}.l; %predshift takes up way too much storage, recompute each time...
end
clear k tmp

% reverse map, contains clusterindex for each coefficient
revmap = repelem(1:numel(kgx), 1, n);

%% main attack loop
parfor i = 1:2; x = 1; end %just to start parallel pool

tattack = tic;

for iter = 1:maxiterations
    fprintf('Iteration: %d/%d\n', iter, maxiterations);
    iterations(sweep_numfaults, sweep_repeat) = iter;
    tic
    
    % update noiserange (needed in case of mergers)
    m = max(cellfun(@(x) numel(x.idx), kgx));
    m = min(m, 600);
    r = m*eta*vals(end);
    noiserange = -r:r;
    
    fprintf('  1: Updating expected error distribution\n');
    ppm = ParforProgMon('Progress ', numel(kgx));
    parfor k = 1:numel(kgx)%can be parallelized
        predshift = kgx{k}.kg*kgx{k}.l;
        kgx{k}.dist = zeros(numel(noiserange), numfaults, 'single');
        for i = 1:numfaults
            subs = predshift(:, i) - noiserange(1) + 1;
            kgx{k}.dist(:, i) = accumarray(subs, kgx{k}.pkg, [numel(noiserange), 1]);
        end
        ppm.increment(); %#ok<PFBNS>
    end
    clear ppm k predshift i subs 
    
    
    
    fprintf('  2: Combining all error distributions\n');

    r = 600;
    nr = -r:r;
    
    
    %maximum number of faults to process at once
    max_junksize = 500;
    numchunks = ceil(numfaults/max_junksize);
    fdist = cell(numel(kgx), numchunks);
    for chunkidx = 1:numchunks
        chunkidxes = ((chunkidx-1)*max_junksize +1):(chunkidx*max_junksize);
        chunkidxes = chunkidxes(chunkidxes <= numfaults);
        
        currnum = numel(chunkidxes);
        
        
        s = struct('up', complex(zeros(2*numel(nr), currnum, 'single')), 'down', complex(zeros(2*numel(nr), currnum, 'single')));

        % prepare a tree data structure
        % needed layers: 1: k, 2: (s,r), then for all combinations, starting from n
        layers = floor(log2(numel(kgx)));
        tree = cell(layers, 1);
        tmp = numel(kgx);
        for i = 1:layers
            t = cell(tmp, 1);
            t(:) = {s};
            tree{i} = t; 
            tmp = ceil(tmp/2);
        end
        clear tmp t
        
        % initialize bottom layer
        rt = (numel(noiserange)-1)/2;
        x = 2^nextpow2(2*numel(nr)) - numel(noiserange);
        for i = 1:numel(kgx)
            t = s;
            t.up = fft([kgx{i}.dist(rt+1:end, chunkidxes); zeros(x, currnum, 'single'); kgx{i}.dist(1:rt, chunkidxes)], [], 1);
            tree{1}{i} = t;
        end
        clear x i t rt
        
        % walk up the tree
        for layer = 2:layers
            for i = 1:numel(tree{layer})
                idx = (i-1)*2 + 1;
                if idx + 1 > numel(tree{layer-1})
                    tree{layer}{i}.up = tree{layer-1}{idx}.up;
                else
                    tree{layer}{i}.up = tree{layer-1}{idx}.up.*tree{layer-1}{idx+1}.up;
                end
            end
        end
        clear layer i idx

        % now walk down the tree again using the dist of all others
        %top layer
        for i = 1:numel(tree{end})
            init = false;
            for j = 1:numel(tree{end})
                if i == j
                    continue
                end
                if ~init
                    tree{end}{i}.down = tree{end}{j}.up;
                    init = true;
                else
                    tree{end}{i}.down = tree{end}{i}.down.*tree{end}{j}.up;
                end
            end
        end
        clear i j init

        %lower layers
        for layer = (layers-1):-1:1
            for i = 1:numel(tree{layer})
                sibling = bitxor((i-1), 1) + 1;
                parent = floor((i - 1)/2) + 1;
                if(sibling > numel(tree{layer}))
                    tree{layer}{i}.down = tree{layer+1}{parent}.down;
                else
                    tree{layer}{i}.down = tree{layer+1}{parent}.down.*tree{layer}{sibling}.up;
                end
            end
        end
        clear layer i sibling parent
        
        for k = 1:numel(kgx)
            cdist = abs(ifft(tree{1}{k}.down,[], 1, 'symmetric')); %abs just for safety
            idx = [(-(r-1):0) + size(cdist, 1), (0:r)+1];
            cdist = cdist(idx, :);
            
            fdist{k, chunkidx} = single(cdist);
        end
        clear k cdist idx
        
        clear tree
    end
    
    clear numchunks chunkidx chunkidxs currnum    
    
    % now do the bayesian updating for all key parts
    fprintf('  3: Update probability of all key candidates\n');
    for k = 1:numel(kgx) %parfor makes this a lot slower somehow...
%         pkgcurr = kgx{k}.pkg; % must not use previous result as priors!
        kgix = kgx{k}.kg - vals(1) + 1;
        pkgcurr = prod(pbino(kgix), 2);
        pkgcurr = log(pkgcurr);
        
        cdist = double(cell2mat(fdist(k, :)));
        
        predshift = kgx{k}.kg*kgx{k}.l;
        
        % new vectorized version
        premcdf = cumsum(cdist, 1);
        
        cutoffidx = cutoff - nr(1) + 1;
        predidx = bsxfun(@minus, cutoffidx, predshift);
        ns2 = size(predshift, 2);
        predidx = bsxfun(@plus, predidx, (0:(ns2-1))*numel(nr));
        predcdf = premcdf(predidx);
        
        % TODO could introduce like an "error probability", i.e., say that we expect that, e.g., 5% of the ineqs are wrong
        %  would boil down to making a weighted sum of left and right (x*left + (1-x)*right)
        predcdf(:, ~left) = 1 - predcdf(:, ~left);
        predcdf = log(predcdf);
        pkgcurr = pkgcurr + sum(predcdf, 2);
        pkgcurr = pkgcurr - max(pkgcurr);
        
        pkgcurr = exp(pkgcurr);
        
        pkgcurr = pkgcurr/sum(pkgcurr);
        kgx{k}.pkg = pkgcurr;
    end
    clear k i cdist dist prem premcdf cutoffidx predidx predcdf pkgcurr kgix ns2 fdist
    
    
    % marginalize
    coeffprobs = zeros(2*kyber, nx);
    for i = 1:numel(kgx)
        pkg = kgx{i}.pkg;
        for j = 1:numel(kgx{i}.idx)
            idx = kgx{i}.idx(j);
            coeffprobs(idx, :) = accumarray(kgx{i}.kg(:, j) - vals(1) + 1, pkg, [numel(vals), 1]);
        end
    end
    coeffentropy = log2(coeffprobs);
    coeffentropy(coeffprobs == 0) = 0;
    coeffentropy = -sum(coeffprobs.*coeffentropy, 2);
    
    [maxprobs, midx] = max(coeffprobs, [], 2);
    recovered = midx + vals(1) - 1;
    correct = recovered' == se;
    
    [~, si] = sort(maxprobs, 'descend');
    firsterr = find(~correct(si), 1, 'first');
    if isempty(firsterr)
        firsterr = 2*kyber + 1; % no errors left --> set firsterr to 1 past maximum
    end
    si = si(1:kyber);
    si = sort(si);
    
    numerr = sum(~(correct(si)));
    
    clear i j si midx recovered
    
    fprintf('  4: Have %d incorrect coefficients\n', numerr);
    fprintf('     First error after %d most likely coeffs\n', firsterr);
    firsterror_best(sweep_numfaults, sweep_repeat) = max(firsterr, firsterror_best(sweep_numfaults, sweep_repeat));
    firsterror_iterations(sweep_numfaults, sweep_repeat, iter) = firsterr;
    % plot stuff
    
    cprob = coeffprobs(sub2ind(size(coeffprobs), 1:2*kyber, se - vals(1) + 1));
    rank = sum(bsxfun(@ge, coeffprobs, cprob'), 2);
    
    if(plots)
        figure('name', sprintf('Iteration %d', iter));
        subplot(4, 1, 1);
        plot(1:2*kyber, rank);
        title('rank of correct key');
        xlim([1, 2*kyber]);
        ylim([1, numel(vals)]);
        subplot(4, 1, 2);
        plot(1:2*kyber, cprob);
        title('prob of correct key');
        ylim([0 1]);
        xlim([1, 2*kyber]);
        subplot(4, 1, 3);
        histogram(rank, (1:numel(vals)+1)-0.5, 'normalization', 'probability');
        title('hist of key rank');
        ylim([0 1]);
    end
    
    if numerr == 0
        fprintf('Attack successful!\n');
        successes(sweep_numfaults, sweep_repeat) = true;
        break;
    end
    
    % determine how many keys we can discard
    if domerge

        pplain = cellfun(@(x) {x.pkg}, kgx);
        numpkg = sum(cellfun(@(x) numel(x), pplain));
        totnumkg = sum(log2(cellfun(@(x) numel(x), pplain)));
        [used, perr] = selectKeys(pplain, err_threshold);
        clear pplain

        if(plots)
            subplot(4, 1, 4);
            plot(perr);
            xlim([1 numel(perr)]);
            title('proberr');
        end



        numexcluded = find(perr < err_threshold, 1, 'last');
        fprintf('     can exclude %d/%d (%.2f%%) of key guesses\n', numexcluded, numpkg, (numexcluded/numpkg*100));

        newnumkg = sum(log2(cellfun(@(x) sum(x), used)));
        fprintf('     reduction: 2^%.0f --> 2^%.0f possible keys\n', totnumkg, newnumkg);

        % compare this to a coefficient-wise reduction
        cc = mat2cell(coeffprobs, ones(size(coeffprobs, 1), 1));
        [ccused, ccperr] = selectKeys(cc, err_threshold);
        cctotnumkg = sum(log2(cellfun(@(x) numel(x), kgxcoeff)));
        ccnewnumkg = sum(log2(cellfun(@(x) sum(x), ccused)));
        ccnumexcluded = find(ccperr < err_threshold, 1, 'last');
        ccnumpkg = sum(cellfun(@(x) numel(x), cc));
        fprintf('     can exclude %d/%d (%.2f%%) of key guesses\n', ccnumexcluded, ccnumpkg, (ccnumexcluded/ccnumpkg*100));
        fprintf('     reduction: 2^%.0f --> 2^%.0f possible keys (coeff wise)\n', cctotnumkg, ccnewnumkg);

        condentropy = zeros(size(coeffentropy));
        mutinfo = zeros(size(coeffentropy));
        for k = 1:numel(mutinfo)
            clusteridx = revmap(k);
            inclusteridx = kgx{clusteridx}.idx == k;

            % compute the entropy of the cluster without the coeff --> marginalize
            [~, pymargin] = marginalizeCluster(kgx{clusteridx}, ~inclusteridx);
            hmargin = entropy(pymargin);
            mutinfo(k) = hmargin - condentropy(k);
        end
        clear k clusteridx inclusteridx kg pkg pc ic pymargin
    end
    
    if domerge
        
        failed = false;
        for k = 1:numel(kgx)
            idx = used{k};
            kgx{k}.pkg = kgx{k}.pkg(idx);
            kgx{k}.pkg = kgx{k}.pkg/sum(kgx{k}.pkg);
            kgx{k}.kg = kgx{k}.kg(idx, :);
            
            if ~ismember(kgx{k}.kg, se(kgx{k}.idx), 'rows')
                failed = true;
                break
            end
        end
        if failed
            fprintf('Correct key discarded, attack failed!\n');
            break;
        end
        
        usedkgx = false(size(kgx));
    
        openvictim = false;
        currvictim = 0;
        nextvictimidx = 0;
        newkgx = cell(0);

        while ~all(usedkgx)
            %get the next unused kgx
            nextidx = find(~usedkgx, 1, 'first');

            nextcluster = kgx{nextidx};
            usedkgx(nextidx) = true;

            numkg = size(nextcluster.kg, 1);
            mergeidx = [];

            while true
                %find a victim cluster to break apart
                if ~openvictim
                    currvictim = find(~usedkgx, 1, 'first');
                    if isempty(currvictim) % ran out of clusters
                        break;
                    end
                    nextvictimidx = 1;
                end

                newkgs = numel(unique(kgx{currvictim}.kg(:, nextvictimidx)));
                numkg = numkg * newkgs;

                if numkg > mergethreshold; break; end % if the new number of key guesses is too large --> abort

                openvictim = true;
                usedkgx(currvictim) = true;
                mergeidx(end+1) = nextvictimidx;

                %check if there are more indizes for the current victim
                if (nextvictimidx + 1) > numel(kgx{currvictim}.idx)
                    openvictim = false;
                    nextcluster = mergeClusters(nextcluster, [], kgx{currvictim}, mergeidx, l);
                    mergeidx = [];
                else
                    nextvictimidx = nextvictimidx + 1;
                end
            end
            if not(isempty(mergeidx))
                nextcluster = mergeClusters(nextcluster, [], kgx{currvictim}, mergeidx, l);
            end

            newkgx{end+1} = nextcluster;
        end

        if openvictim
            nextcluster = kgx{currvictim};
            tmp = nextvictimidx:numel(nextcluster.idx);
            nextcluster.idx = nextcluster.idx(tmp);
            [nextcluster.kg, nextcluster.pkg] = marginalizeCluster(nextcluster, tmp);
            nextcluster.l =  l(nextcluster.idx, :);
            newkgx{end+1} = nextcluster;
        end
        kgx = newkgx;
        clear usedkgx openvictim currvictim nextvictimidx nextidx nextcluster numkg mergeidx tmp newkgx
        
        %update revmap
        for k = 1:numel(kgx)
            revmap(kgx{k}.idx) = k;
        end
    end
    
    toc
end

fprintf('Attack time: \n');
toc(tattack)

%%

    end
    
    if store
        save result_newhope512.mat successes iterations firsterror_best firsterror_iterations repeats faultrange 
    end
end

%% analysis

if numel(faultrange) > 1
    sr = mean(successes, 2);
    figure;
    plot(faultrange, sr);
end


