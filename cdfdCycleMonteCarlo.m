function MC = cdfdCycleMonteCarlo(C, varargin)
% cdfdCycleMonteCarlo
% Streaming Monte Carlo estimator for probabilistic CDFD cycle coefficients.
%
% PURPOSE
%   Estimate the Qian--Qian / Kalpazidou probabilistic cycle representation
%   of a circular component C by simulating the associated Markov chain and
%   applying chronological loop erasure.
%
%   For each strongly connected component a of C, with edge-flow volume S_a,
%   the estimator uses
%
%       q_hat_gamma      = N_gamma / T_a
%       lambda_hat_gamma = S_a * q_hat_gamma
%
%   where N_gamma is the number of times cycle gamma is sampled and T_a is
%   the number of Markov transitions simulated in component a.
%
% INPUT
%   C   n x n non-negative circular component, typically C_BFF.
%
% OPTIONS
%   'Tol'              default 1e-12
%   'BudgetMode'       default 'cycles'. One of:
%                        'cycles'       simulate NumSamples cycles per component
%                        'transitions'  simulate at least NumTransitions
%                                       Markov transitions per component
%   'NumSamples'       default 10000. Used when BudgetMode='cycles'.
%   'NumTransitions'   default 100000. Used when BudgetMode='transitions'.
%   'Seed'             default []. If non-empty, calls rng(Seed).
%   'TopK'             default Inf. Number of ranked observed cycles to return.
%   'RankBy'           default 'throughput'. Passed to cdfdTopCycles.
%   'CheckBalance'     default true.
%   'Verbose'          default true.
%   'MaxTransitions'   default Inf. Safety cap per component.
%
% OUTPUT
%   MC   struct with fields:
%
%        CycleTable       estimated observed-cycle table with columns:
%                           Component
%                           Cycle
%                           Length
%                           Count
%                           q
%                           lambda
%                           S
%                           NumTransitions
%                           NumSamples
%                           ComponentVolumeShare
%
%        TopCycles        ranked observed cycles using cdfdTopCycles
%
%        Estimates        struct with aggregate Monte Carlo estimates:
%                           V_C
%                           T_C
%                           Lbar_C
%                           VolumeFromCycles
%                           VolumeGap
%                           VolumeRelError
%                           NumSamples
%                           NumTransitions
%                           NumObservedCycles
%
%        ComponentTable   per-component simulation diagnostics
%        Options          parsed options
%
% NOTES
%   - This is a Monte Carlo estimator, not an exact all-cycle enumerator.
%   - CycleTable contains only observed cycles.
%   - Passing MC.CycleTable to cdfdCycleSummaries will summarise the observed
%     estimated table, not a guaranteed complete exact decomposition.
%   - Self-loops are allowed and are sampled as length-one cycles.

    opts = parseMonteCarloOptions(varargin{:});
    tol = opts.Tol;

    validateSquareNonnegative(C, tol);

    if ~isempty(opts.Seed)
        rng(opts.Seed);
    end

    C = sparse(double(C));
    C(abs(C) <= tol) = 0;
    C = sparse(C);

    n = size(C, 1);
    V_C = full(sum(C(:)));

    if opts.CheckBalance
        balanceResidual = full(max(abs(sum(C, 2) - sum(C, 1).')));
        scale = full(max(1, max(sum(C, 2))));

        if balanceResidual > 100 * tol * scale
            warning('C does not appear balanced. Max balance residual = %.3g.', balanceResidual);
        end
    end

    components = getPositiveSupportComponents(C, tol);

    if isempty(components)
        MC = emptyMonteCarloOutput(n, V_C, opts);
        return;
    end

    rowsComponent = zeros(0, 1);
    rowsCycle = cell(0, 1);
    rowsLength = zeros(0, 1);
    rowsCount = zeros(0, 1);
    rowsQ = zeros(0, 1);
    rowsLambda = zeros(0, 1);
    rowsS = zeros(0, 1);
    rowsTransitions = zeros(0, 1);
    rowsSamples = zeros(0, 1);
    rowsComponentVolumeShare = zeros(0, 1);

    compRows = {};

    totalSamples = 0;
    totalTransitions = 0;

    if opts.Verbose
        fprintf('\nRunning CDFD cycle Monte Carlo on %d component(s)...\n', numel(components));
    end

    for a = 1:numel(components)
        nodes = components{a};
        Csub = C(nodes, nodes);

        componentResult = simulateComponent(Csub, nodes, a, V_C, opts);

        totalSamples = totalSamples + componentResult.NumSamples;
        totalTransitions = totalTransitions + componentResult.NumTransitions;

        Tlocal = componentResult.CycleTable;

        if ~isempty(Tlocal)
            rowsComponent = [rowsComponent; Tlocal.Component]; %#ok<AGROW>
            rowsCycle = [rowsCycle; Tlocal.Cycle]; %#ok<AGROW>
            rowsLength = [rowsLength; Tlocal.Length]; %#ok<AGROW>
            rowsCount = [rowsCount; Tlocal.Count]; %#ok<AGROW>
            rowsQ = [rowsQ; Tlocal.q]; %#ok<AGROW>
            rowsLambda = [rowsLambda; Tlocal.lambda]; %#ok<AGROW>
            rowsS = [rowsS; Tlocal.S]; %#ok<AGROW>
            rowsTransitions = [rowsTransitions; Tlocal.NumTransitions]; %#ok<AGROW>
            rowsSamples = [rowsSamples; Tlocal.NumSamples]; %#ok<AGROW>
            rowsComponentVolumeShare = [rowsComponentVolumeShare; Tlocal.ComponentVolumeShare]; %#ok<AGROW>
        end

        compRows(end+1, :) = { ...
            a, ...
            numel(nodes), ...
            nnz(Csub > tol), ...
            componentResult.S, ...
            componentResult.ComponentVolumeShare, ...
            componentResult.NumSamples, ...
            componentResult.NumTransitions, ...
            componentResult.NumObservedCycles, ...
            componentResult.T_C, ...
            componentResult.Lbar_C, ...
            componentResult.VolumeFromCycles, ...
            componentResult.VolumeRelError ...
        }; %#ok<AGROW>

        if opts.Verbose
            fprintf(['Component %d: nodes=%d, S=%.6g, samples=%d, ' ...
                     'transitions=%d, observed cycles=%d, T_C_hat=%.6g\n'], ...
                a, ...
                numel(nodes), ...
                componentResult.S, ...
                componentResult.NumSamples, ...
                componentResult.NumTransitions, ...
                componentResult.NumObservedCycles, ...
                componentResult.T_C);
        end
    end

    CycleTable = table( ...
        rowsComponent, ...
        rowsCycle, ...
        rowsLength, ...
        rowsCount, ...
        rowsQ, ...
        rowsLambda, ...
        rowsS, ...
        rowsTransitions, ...
        rowsSamples, ...
        rowsComponentVolumeShare, ...
        'VariableNames', { ...
            'Component', ...
            'Cycle', ...
            'Length', ...
            'Count', ...
            'q', ...
            'lambda', ...
            'S', ...
            'NumTransitions', ...
            'NumSamples', ...
            'ComponentVolumeShare' ...
        } ...
    );

    if ~isempty(CycleTable)
        CycleTable = sortrows(CycleTable, {'lambda', 'Length'}, {'descend', 'ascend'});
    end

    if isempty(compRows)
        ComponentTable = table();
    else
        ComponentTable = cell2table(compRows, ...
            'VariableNames', { ...
                'Component', ...
                'NumNodes', ...
                'NumEdges', ...
                'S', ...
                'ComponentVolumeShare', ...
                'NumSamples', ...
                'NumTransitions', ...
                'NumObservedCycles', ...
                'T_C', ...
                'Lbar_C', ...
                'VolumeFromCycles', ...
                'VolumeRelError' ...
            });
    end

    T_C_hat = full(sum(CycleTable.lambda));
    volumeFromCycles = full(sum(CycleTable.lambda .* CycleTable.Length));
    volumeGap = volumeFromCycles - V_C;
    volumeRelError = abs(volumeGap) / max(1, abs(V_C));

    if T_C_hat > tol
        Lbar_hat = V_C / T_C_hat;
    else
        Lbar_hat = NaN;
    end

    if isempty(CycleTable)
        TopCycles = table();
    else
        TopCycles = cdfdTopCycles(CycleTable, ...
            'TopK', opts.TopK, ...
            'RankBy', opts.RankBy, ...
            'Tol', tol);
    end

    Estimates = struct();
    Estimates.V_C = V_C;
    Estimates.T_C = T_C_hat;
    Estimates.Lbar_C = Lbar_hat;
    Estimates.VolumeFromCycles = volumeFromCycles;
    Estimates.VolumeGap = volumeGap;
    Estimates.VolumeRelError = volumeRelError;
    Estimates.NumSamples = totalSamples;
    Estimates.NumTransitions = totalTransitions;
    Estimates.NumObservedCycles = height(CycleTable);
    Estimates.NumComponents = numel(components);

    MC = struct();
    MC.CycleTable = CycleTable;
    MC.TopCycles = TopCycles;
    MC.Estimates = Estimates;
    MC.ComponentTable = ComponentTable;
    MC.Options = opts;
end


function result = simulateComponent(Csub, globalNodes, componentIndex, V_C, opts)
% simulateComponent
% Run loop-erased cycle sampling inside one strongly connected component.

    tol = opts.Tol;

    Csub = full(Csub);
    Csub(abs(Csub) <= tol) = 0;

    m = size(Csub, 1);

    s = full(sum(Csub, 2));
    S = full(sum(s));

    if S <= tol
        error('Component %d has zero circular edge-flow volume.', componentIndex);
    end

    if any(s <= tol)
        error('Component %d contains a node with zero outflow.', componentIndex);
    end

    P = Csub ./ s;
    startProb = s / S;

    componentVolumeShare = S / max(V_C, eps);

    path = sampleDiscrete(startProb);

    cycleCount = containers.Map('KeyType', 'char', 'ValueType', 'double');
    cycleMap = containers.Map('KeyType', 'char', 'ValueType', 'any');

    numSamples = 0;
    numTransitions = 0;

    budgetMode = lower(string(opts.BudgetMode));

    while true
        switch budgetMode
            case "cycles"
                if numSamples >= opts.NumSamples
                    break;
                end

            case "transitions"
                if numTransitions >= opts.NumTransitions
                    break;
                end

            otherwise
                error('BudgetMode must be cycles or transitions.');
        end

        if numTransitions >= opts.MaxTransitions
            error('Component %d reached MaxTransitions=%g before completing requested budget.', ...
                componentIndex, opts.MaxTransitions);
        end

        current = path(end);
        next = sampleNextState(P, current, tol);

        numTransitions = numTransitions + 1;

        previousPosition = find(path == next, 1, 'first');

        if isempty(previousPosition)
            path(end+1) = next; %#ok<AGROW>
        else
            cycleLocal = path(previousPosition:end);
            path = path(1:previousPosition);

            cycleGlobal = globalNodes(cycleLocal);
            cycleCanonical = canonicalDirectedCycle(cycleGlobal);
            key = cycleKey(cycleCanonical);

            if isKey(cycleCount, key)
                cycleCount(key) = cycleCount(key) + 1;
            else
                cycleCount(key) = 1;
                cycleMap(key) = cycleCanonical;
            end

            numSamples = numSamples + 1;
        end
    end

    keys = cycleCount.keys;
    numObserved = numel(keys);

    rowsComponent = zeros(numObserved, 1);
    rowsCycle = cell(numObserved, 1);
    rowsLength = zeros(numObserved, 1);
    rowsCount = zeros(numObserved, 1);
    rowsQ = zeros(numObserved, 1);
    rowsLambda = zeros(numObserved, 1);
    rowsS = zeros(numObserved, 1);
    rowsTransitions = zeros(numObserved, 1);
    rowsSamples = zeros(numObserved, 1);
    rowsComponentVolumeShare = zeros(numObserved, 1);

    for r = 1:numObserved
        key = keys{r};
        cyc = cycleMap(key);
        count = cycleCount(key);

        qhat = count / numTransitions;
        lambdaHat = S * qhat;

        rowsComponent(r) = componentIndex;
        rowsCycle{r} = cyc;
        rowsLength(r) = numel(cyc);
        rowsCount(r) = count;
        rowsQ(r) = qhat;
        rowsLambda(r) = lambdaHat;
        rowsS(r) = S;
        rowsTransitions(r) = numTransitions;
        rowsSamples(r) = numSamples;
        rowsComponentVolumeShare(r) = componentVolumeShare;
    end

    Tlocal = table( ...
        rowsComponent, ...
        rowsCycle, ...
        rowsLength, ...
        rowsCount, ...
        rowsQ, ...
        rowsLambda, ...
        rowsS, ...
        rowsTransitions, ...
        rowsSamples, ...
        rowsComponentVolumeShare, ...
        'VariableNames', { ...
            'Component', ...
            'Cycle', ...
            'Length', ...
            'Count', ...
            'q', ...
            'lambda', ...
            'S', ...
            'NumTransitions', ...
            'NumSamples', ...
            'ComponentVolumeShare' ...
        } ...
    );

    if ~isempty(Tlocal)
        Tlocal = sortrows(Tlocal, {'lambda', 'Length'}, {'descend', 'ascend'});
    end

    T_C_hat = full(sum(Tlocal.lambda));
    volumeFromCycles = full(sum(Tlocal.lambda .* Tlocal.Length));

    if T_C_hat > tol
        Lbar_hat = S / T_C_hat;
    else
        Lbar_hat = NaN;
    end

    volumeRelError = abs(volumeFromCycles - S) / max(1, abs(S));

    result = struct();
    result.CycleTable = Tlocal;
    result.S = S;
    result.ComponentVolumeShare = componentVolumeShare;
    result.NumSamples = numSamples;
    result.NumTransitions = numTransitions;
    result.NumObservedCycles = numObserved;
    result.T_C = T_C_hat;
    result.Lbar_C = Lbar_hat;
    result.VolumeFromCycles = volumeFromCycles;
    result.VolumeRelError = volumeRelError;
end


function components = getPositiveSupportComponents(C, tol)
% getPositiveSupportComponents
% Return components of the positive support with positive edge-flow volume.

    n = size(C, 1);

    A = C > tol;

    if nnz(A) == 0
        components = {};
        return;
    end

    G = digraph(A);
    comp = conncomp(G, 'Type', 'strong');
    compIds = unique(comp);

    components = {};

    for k = 1:numel(compIds)
        nodes = find(comp == compIds(k));
        Csub = C(nodes, nodes);
        S = full(sum(Csub(:)));

        if S > tol
            components{end+1, 1} = nodes; %#ok<AGROW>
        end
    end
end


function idx = sampleDiscrete(prob)

    prob = full(prob(:));
    prob(prob < 0) = 0;

    s = sum(prob);

    if s <= 0
        error('Cannot sample from a zero probability vector.');
    end

    prob = prob / s;

    u = rand();
    cs = cumsum(prob);

    idx = find(u <= cs, 1, 'first');

    if isempty(idx)
        idx = find(prob > 0, 1, 'last');
    end
end


function next = sampleNextState(P, current, tol)

    probs = full(P(current, :));
    probs(abs(probs) <= tol) = 0;

    s = sum(probs);

    if s <= tol
        error('Current node %d has zero outgoing transition probability.', current);
    end

    probs = probs / s;

    u = rand();
    cs = cumsum(probs);

    next = find(u <= cs, 1, 'first');

    if isempty(next)
        next = find(probs > 0, 1, 'last');
    end
end


function canonical = canonicalDirectedCycle(cyc)
% canonicalDirectedCycle
% Canonicalise a directed cycle up to cyclic rotation, preserving orientation.

    cyc = cyc(:).';
    L = numel(cyc);

    if L <= 1
        canonical = cyc;
        return;
    end

    canonical = cyc;

    for shift = 1:(L-1)
        rotated = [cyc((shift+1):end) cyc(1:shift)];

        if lexLess(rotated, canonical)
            canonical = rotated;
        end
    end
end


function tf = lexLess(a, b)

    tf = false;

    for k = 1:numel(a)
        if a(k) < b(k)
            tf = true;
            return;
        elseif a(k) > b(k)
            tf = false;
            return;
        end
    end
end


function key = cycleKey(cyc)

    key = strjoin(string(cyc), "_");
end


function MC = emptyMonteCarloOutput(n, V_C, opts)

    CycleTable = table();

    Estimates = struct();
    Estimates.V_C = V_C;
    Estimates.T_C = NaN;
    Estimates.Lbar_C = NaN;
    Estimates.VolumeFromCycles = NaN;
    Estimates.VolumeGap = NaN;
    Estimates.VolumeRelError = NaN;
    Estimates.NumSamples = 0;
    Estimates.NumTransitions = 0;
    Estimates.NumObservedCycles = 0;
    Estimates.NumComponents = 0;

    MC = struct();
    MC.CycleTable = CycleTable;
    MC.TopCycles = table();
    MC.Estimates = Estimates;
    MC.ComponentTable = table();
    MC.Options = opts;
    MC.NumNodes = n;
end


function validateSquareNonnegative(C, tol)

    if ~ismatrix(C) || size(C, 1) ~= size(C, 2)
        error('C must be a square matrix.');
    end

    if any(~isfinite(C(:)))
        error('C contains non-finite entries.');
    end

    if any(C(:) < -tol)
        error('C must be non-negative up to numerical tolerance.');
    end
end


function opts = parseMonteCarloOptions(varargin)

    opts = struct();
    opts.Tol = 1e-12;
    opts.BudgetMode = 'cycles';
    opts.NumSamples = 10000;
    opts.NumTransitions = 100000;
    opts.Seed = [];
    opts.TopK = Inf;
    opts.RankBy = 'throughput';
    opts.CheckBalance = true;
    opts.Verbose = true;
    opts.MaxTransitions = Inf;

    if mod(numel(varargin), 2) ~= 0
        error('Options must be name-value pairs.');
    end

    for k = 1:2:numel(varargin)
        name = varargin{k};
        value = varargin{k+1};

        switch lower(name)
            case 'tol'
                opts.Tol = value;

            case 'budgetmode'
                opts.BudgetMode = value;

            case 'numsamples'
                opts.NumSamples = value;

            case 'numtransitions'
                opts.NumTransitions = value;

            case 'seed'
                opts.Seed = value;

            case 'topk'
                opts.TopK = value;

            case 'rankby'
                opts.RankBy = value;

            case 'checkbalance'
                opts.CheckBalance = value;

            case 'verbose'
                opts.Verbose = value;

            case 'maxtransitions'
                opts.MaxTransitions = value;

            otherwise
                error('Unknown option: %s.', name);
        end
    end

    validatePositiveIntegerLike(opts.NumSamples, 'NumSamples');
    validatePositiveIntegerLike(opts.NumTransitions, 'NumTransitions');

    if ~isinf(opts.MaxTransitions)
        validatePositiveIntegerLike(opts.MaxTransitions, 'MaxTransitions');
    end
end


function validatePositiveIntegerLike(x, name)

    if ~isnumeric(x) || ~isscalar(x) || ~isfinite(x) || x < 1 || x ~= round(x)
        error('%s must be a positive integer scalar.', name);
    end
end