%% choose parameter set and attack settings

 %#ok<*UNRCH>

kyberk = 2;
kyber = kyberk*256;

n = 4;

faultrange = 4000:500:10000; %k = 2
% faultrange = 6000:500:12000; % k = 3 (no need to simulate stuff where k=2 already failed)
% faultrange = 8000:500:15000; % k = 4 (no need to simulate stuff where k=3 already failed)
faultrange = 7500; %for demonstration
numfaults = 10000;
maxiterations = 10;

repeats = 1;

%exclusion and filtering parameters
err_threshold = 0.005;
domerge = false;
mergethreshold = 1000; % maximum cluster size

errate = 0; %allows to introduce an additional error in the (simulated) data

realdata = true;
device = false;
devicetest = 4;

plots = true;
store = false;

%%


if device
    repeats = 1;
    faultrange = 10000;
    store = false;
    plots = true;
    errate = 0; %we don't want to introduce an additional error in real data
    
    fprintf('Real device: test %d\n', devicetest);
end

successes = false(numel(faultrange), repeats);
iterations = zeros(size(successes));
firsterror_best = zeros(size(successes));
firsterror_iterations = zeros(numel(faultrange), repeats, maxiterations);
sweep_repeat = 1;
sweep_numfaults = 1;

%% compute constants
consts %initializes err_round_u and err_round_v

eta = 2;
vals = -eta:eta;

valsx = vals - vals(1);
pbino = binopdf(valsx, 2*eta, 0.5);

binoentropy = -sum(pbino.*log2(pbino));

%distribution of (e1 + du), obtained via convolution
valssum = (err_round_u(1, 1)-eta):(err_round_u(end, 1)+eta);
psum = conv(pbino, err_round_u(:, 2), 'full');

psumcdf = cumsum(psum);


%% sweep start
for sweep_repeat = 1:repeats
    t_repeat = tic;
    for sweep_numfaults = 1:numel(faultrange)
        if ~device
            numfaults = faultrange(sweep_numfaults);
            fprintf('SWEEP. NUMFAULTS: %d, REPEAT: %d\n', numfaults, sweep_repeat)
        end
    

%% generate or read in test data
nx = numel(vals);

if ~realdata && ~device %simulation
    
    s = binornd(nx-1, 0.5, 1, kyber) - eta;
    e = binornd(nx-1, 0.5, 1, kyber) - eta;
    se = [s, e];

    % s gets multiplied with samples from (e1+du) = psum
    l1 = inversionsample(valssum, psumcdf, kyber, numfaults);
    % e gets multiplied with samples from e2 = pbino
    l2 = binornd(nx-1, 0.5, kyber, numfaults) - eta;
    l = [l1 ; l2];
    % the full noise vector
    realnoise = se*l;

    % generate random cutoffs and check for every single one if we are on the left
    coff = 10;
    r = -coff:coff;
    cutoff = r(randi(numel(r), 1, numfaults));
    % cutoff = zeros(1, numfaults);
    % caution: cdf is defined as P(X <= x)
    % but our fault gives us the information whether X < x or X >= x
    % solution: subtract 1 from x, get X <= (x-1) or X > (x-1)
    % for the simulation its not important, but it is when using real data
    left = realnoise <= cutoff; %use <=, so that we are consistent with definition of cdf
    
    clear coff r
else
    
    if device
    %     % ineqs from device
        f1 = sprintf('%d/ineqs_device.txt', devicetest);
        f2 = sprintf('%d/fixed_es.txt', devicetest); 
    else
        
        % kyber, simulated using reference implementation
        f1 = 'ineqs.txt';
        f2 = 'ineqs_es.txt';

        % simulate
        sstring = sprintf('./kyber%d_fault s -n %d', kyber, numfaults);
        system(sstring);

        movefile(f1, ['test_data/', f1]);
        movefile(f2, ['test_data/', f2]);

        % masked decoder, simulated
    %     f1 = 'ineqs_masked.txt';
    %     f2 = 'ineqs_es_masked.txt';

        % ineqs using the fixed key
    %     f1 = 'ineqs_fixed.txt';
    %     f2 = 'fixed_es.txt';
    end
    
    d = readtable(['test_data/', f1], 'FileType', 'text', 'ReadVariableNames', false, 'MultipleDelimsAsOne', true, 'HeaderLines', 0);
    l = table2array(d(:, 1:end-2));
    left = strcmp(table2cell(d(:, end-1)), '<')';
    cutoff = table2array(d(:, end))';
    cutoff = cutoff - 1; % correction mentioned above
    
    numvars = size(l, 2);
    fid = fopen(['test_data/', f2]);
    es = textscan(fid, '%f', numvars);
    es = es{1}';
    fclose(fid);
    
    numfaults = size(l, 1);
    
    %introduce an error
    if errate > 0
        selerr = false(size(left));
        numerr = round(numel(selerr)*errate);
        selerr(1:numerr) = true;
        selerr = selerr(randperm(numel(selerr)));
        %x = binornd(1, errate, size(left)) > 0;
        left(selerr) = ~left(selerr);
        clear numerr
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
end


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

r = n*eta*valssum(end);
noiserange = -r:r;
tmp.dist = zeros(numel(noiserange), numfaults);

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
    m = min(m, 250);
    r = m*eta*valssum(end);
    noiserange = -r:r;
    
    fprintf('  1: Updating expected error distribution\n');
    parfor k = 1:numel(kgx)%can be parallelized
        predshift = kgx{k}.kg*kgx{k}.l;
        kgx{k}.dist = zeros(numel(noiserange), numfaults);
        for i = 1:numfaults
            subs = predshift(:, i) - noiserange(1) + 1;
            kgx{k}.dist(:, i) = accumarray(subs, kgx{k}.pkg, [numel(noiserange), 1]);
        end
    end
    clear k predshift i subs 
    
    
    
    fprintf('  2: Combining all error distributions\n');

    r = 250;
    nr = -r:r;
    
    
    %maximum number of faults to process at once
    max_junksize = 1000;
    numchunks = ceil(numfaults/max_junksize);
    fdist = cell(numel(kgx), numchunks);
    for chunkidx = 1:numchunks
        chunkidxes = ((chunkidx-1)*max_junksize +1):(chunkidx*max_junksize);
        chunkidxes = chunkidxes(chunkidxes <= numfaults);
        
        currnum = numel(chunkidxes);
        
        
        s = struct('up', complex(zeros(2*numel(nr), currnum)), 'down', complex(zeros(2*numel(nr), currnum)));

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
            t.up = fft([kgx{i}.dist(rt+1:end, chunkidxes); zeros(x, currnum); kgx{i}.dist(1:rt, chunkidxes)], [], 1);
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
            
            fdist{k, chunkidx} = cdist;
        end
        clear k cdist idx
        
        clear tree
    end
    
    clear numchunks chunkidx chunkidxs currnum   
    
    % now do the bayesian updating for all key parts
    fprintf('  3: Update probability of all key candidates\n');
    for k = 1:numel(kgx)
        kgix = kgx{k}.kg - vals(1) + 1;
        pkgcurr = prod(pbino(kgix), 2);
        pkgcurr = log(pkgcurr);

        cdist = cell2mat(fdist(k, :));
        
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
        predcdf(predcdf < 0) = 0; %safety measure
        predcdf = log(predcdf);
        sel = max(predcdf, [], 1) > -Inf; %exclude equations where none of the key guesses matches --> must be a wrong equation or something
        pkgcurr = pkgcurr + sum(predcdf(:, sel), 2);
        m = max(pkgcurr(~isinf(pkgcurr)));
        if isempty(m)
            continue;
        end
        pkgcurr = pkgcurr - m; %need to max at 0, otherwise we could have inf - inf = NaN
        
        pkgcurr = exp(pkgcurr);
        
        s = sum(pkgcurr);
        if(s ~= 0) %safety net
            pkgcurr = pkgcurr/sum(pkgcurr);
            kgx{k}.pkg = pkgcurr;
        else
            fprintf('all kgs have zero probability, do not update!\n');
        end
    end
    clear k i cdist dist prem premcdf cutoffidx predidx predcdf pkgcurr kgix ns2 fdist sel
    
    clear tree
    
    clear asdf
    
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
		drawnow;
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
			drawnow;
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
    
    if ~device
        fprintf('Time for one repitition:\n');
        toc(t_repeat)
    end
    if store
        save result.mat successes iterations firsterror_best firsterror_iterations repeats faultrange errate
    end
end

%% analysis

if numel(faultrange) > 1
    sr = mean(successes, 2);
    figure;
    plot(faultrange, sr);
end


