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
%   'BatchTransitions' default []. If non-empty, records fixed-batch
%                      convergence diagnostics every BatchTransitions
%                      transitions. Currently requires BudgetMode='transitions'.
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
%        Estimates        aggregate Monte Carlo estimates:
%                           V_C
%                           T_C
%                           Lbar_C
%                           VolumeFromCycles
%                           VolumeGap
%                           VolumeRelError
%                           NumSamples
%                           NumTransitions
%                           NumObservedCycles
%                           NumComponents
%
%        ComponentTable   per-component simulation diagnostics
%
%        Convergence      convergence diagnostics:
%                           ComponentTrace
%                           Trace
%                           Diagnostics
%
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
    traceTables = {};

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

        if height(Tlocal) > 0
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

        if height(componentResult.Trace) > 0
            traceTables{end+1, 1} = componentResult.Trace; %#ok<AGROW>
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

    if height(CycleTable) > 0
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

    if height(CycleTable) == 0
        TopCycles = table();
    else
        TopCycles = cdfdTopCycles(CycleTable, ...
            'TopK', opts.TopK, ...
            'RankBy', opts.RankBy, ...
            'Tol', tol);
    end

    if isempty(traceTables)
        ComponentTrace = table();
    else
        ComponentTrace = vertcat(traceTables{:});
    end

    [Trace, convergenceDiagnostics] = aggregateComponentTraces( ...
        ComponentTrace, V_C, numel(components), tol);

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

    Convergence = struct();
    Convergence.ComponentTrace = ComponentTrace;
    Convergence.Trace = Trace;
    Convergence.Diagnostics = convergenceDiagnostics;

    MC = struct();
    MC.CycleTable = CycleTable;
    MC.TopCycles = TopCycles;
    MC.Estimates = Estimates;
    MC.ComponentTable = ComponentTable;
    MC.Convergence = Convergence;
    MC.Options = opts;
end


function result = simulateComponent(Csub, globalNodes, componentIndex, V_C, opts)
% simulateComponent
% Run loop-erased cycle sampling inside one strongly connected component.

    tol = opts.Tol;

    Csub = full(Csub);
    Csub(abs(Csub) <= tol) = 0;

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
    totalCycleLengthSum = 0;

    traceEnabled = ~isempty(opts.BatchTransitions);
    batchNumber = 0;
    nextBatchTransition = opts.BatchTransitions;

    batchStartTransitions = 0;
    batchStartSamples = 0;
    batchStartCycleLengthSum = 0;

    traceRows = {};

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
            totalCycleLengthSum = totalCycleLengthSum + numel(cycleCanonical);
        end

        if traceEnabled && numTransitions >= nextBatchTransition
            batchNumber = batchNumber + 1;

            traceRows(end+1, :) = makeComponentTraceRow( ...
                componentIndex, ...
                batchNumber, ...
                S, ...
                componentVolumeShare, ...
                numSamples, ...
                numTransitions, ...
                totalCycleLengthSum, ...
                cycleCount, ...
                batchStartSamples, ...
                batchStartTransitions, ...
                batchStartCycleLengthSum, ...
                tol); %#ok<AGROW>

            batchStartTransitions = numTransitions;
            batchStartSamples = numSamples;
            batchStartCycleLengthSum = totalCycleLengthSum;
            nextBatchTransition = nextBatchTransition + opts.BatchTransitions;
        end
    end

    if traceEnabled
        if batchNumber == 0 || batchStartTransitions < numTransitions
            batchNumber = batchNumber + 1;

            traceRows(end+1, :) = makeComponentTraceRow( ...
                componentIndex, ...
                batchNumber, ...
                S, ...
                componentVolumeShare, ...
                numSamples, ...
                numTransitions, ...
                totalCycleLengthSum, ...
                cycleCount, ...
                batchStartSamples, ...
                batchStartTransitions, ...
                batchStartCycleLengthSum, ...
                tol); %#ok<AGROW>
        end
    end

    Tlocal = makeCycleTableFromCounts( ...
        cycleCount, ...
        cycleMap, ...
        componentIndex, ...
        S, ...
        numTransitions, ...
        numSamples, ...
        componentVolumeShare);

    T_C_hat = full(sum(Tlocal.lambda));
    volumeFromCycles = full(sum(Tlocal.lambda .* Tlocal.Length));

    if T_C_hat > tol
        Lbar_hat = S / T_C_hat;
    else
        Lbar_hat = NaN;
    end

    volumeRelError = abs(volumeFromCycles - S) / max(1, abs(S));

    if isempty(traceRows)
        Trace = table();
    else
        Trace = cell2table(traceRows, ...
            'VariableNames', { ...
                'Component', ...
                'Batch', ...
                'S', ...
                'ComponentVolumeShare', ...
                'CumulativeSamples', ...
                'CumulativeTransitions', ...
                'ObservedCycles', ...
                'T_C', ...
                'Lbar_C', ...
                'VolumeFromCycles', ...
                'VolumeRelError', ...
                'BatchSamples', ...
                'BatchTransitions', ...
                'Batch_T_C', ...
                'Batch_Lbar_C', ...
                'BatchVolumeFromCycles', ...
                'BatchVolumeRelError' ...
            });
    end

    result = struct();
    result.CycleTable = Tlocal;
    result.Trace = Trace;
    result.S = S;
    result.ComponentVolumeShare = componentVolumeShare;
    result.NumSamples = numSamples;
    result.NumTransitions = numTransitions;
    result.NumObservedCycles = height(Tlocal);
    result.T_C = T_C_hat;
    result.Lbar_C = Lbar_hat;
    result.VolumeFromCycles = volumeFromCycles;
    result.VolumeRelError = volumeRelError;
end


function row = makeComponentTraceRow(componentIndex, batchNumber, S, componentVolumeShare, ...
    numSamples, numTransitions, totalCycleLengthSum, cycleCount, ...
    batchStartSamples, batchStartTransitions, batchStartCycleLengthSum, tol)
% makeComponentTraceRow
% Build one per-component convergence trace row.

    observedCycles = cycleCount.Count;

    if numTransitions > 0
        T_C = S * numSamples / numTransitions;
        volumeFromCycles = S * totalCycleLengthSum / numTransitions;
    else
        T_C = NaN;
        volumeFromCycles = NaN;
    end

    if T_C > tol
        Lbar_C = S / T_C;
    else
        Lbar_C = NaN;
    end

    volumeRelError = abs(volumeFromCycles - S) / max(1, abs(S));

    batchSamples = numSamples - batchStartSamples;
    batchTransitions = numTransitions - batchStartTransitions;
    batchCycleLengthSum = totalCycleLengthSum - batchStartCycleLengthSum;

    if batchTransitions > 0
        batchT_C = S * batchSamples / batchTransitions;
        batchVolumeFromCycles = S * batchCycleLengthSum / batchTransitions;
    else
        batchT_C = NaN;
        batchVolumeFromCycles = NaN;
    end

    if batchT_C > tol
        batchLbar_C = S / batchT_C;
    else
        batchLbar_C = NaN;
    end

    batchVolumeRelError = abs(batchVolumeFromCycles - S) / max(1, abs(S));

    row = { ...
        componentIndex, ...
        batchNumber, ...
        S, ...
        componentVolumeShare, ...
        numSamples, ...
        numTransitions, ...
        observedCycles, ...
        T_C, ...
        Lbar_C, ...
        volumeFromCycles, ...
        volumeRelError, ...
        batchSamples, ...
        batchTransitions, ...
        batchT_C, ...
        batchLbar_C, ...
        batchVolumeFromCycles, ...
        batchVolumeRelError ...
    };
end


function Tlocal = makeCycleTableFromCounts(cycleCount, cycleMap, componentIndex, S, ...
    numTransitions, numSamples, componentVolumeShare)
% makeCycleTableFromCounts
% Convert cycle-count maps into an estimated cycle table.

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

    if height(Tlocal) > 0
        Tlocal = sortrows(Tlocal, {'lambda', 'Length'}, {'descend', 'ascend'});
    end
end


function [Trace, diagnostics] = aggregateComponentTraces(ComponentTrace, V_C, numComponents, tol)
% aggregateComponentTraces
% Aggregate component-level trace rows into network-level batch trace.

    diagnostics = struct();
    diagnostics.Enabled = height(ComponentTrace) > 0;
    diagnostics.NumComponentTraceRows = height(ComponentTrace);
    diagnostics.NumBatches = 0;
    diagnostics.FinalRelChangeT_C = NaN;
    diagnostics.FinalRelChangeLbar_C = NaN;
    diagnostics.MCSE_T_C = NaN;
    diagnostics.MCSE_Lbar_C = NaN;
    diagnostics.StopReason = "fixed_budget";

    if height(ComponentTrace) == 0
        Trace = table();
        return;
    end

    batches = unique(ComponentTrace.Batch);
    rows = {};

    for k = 1:numel(batches)
        b = batches(k);
        sub = ComponentTrace(ComponentTrace.Batch == b, :);

        % Only aggregate batches for which all components have a row.
        if height(sub) ~= numComponents
            continue;
        end

        cumulativeSamples = full(sum(sub.CumulativeSamples));
        cumulativeTransitions = full(sum(sub.CumulativeTransitions));
        observedCycles = full(sum(sub.ObservedCycles));

        T_C = full(sum(sub.T_C));

        if T_C > tol
            Lbar_C = V_C / T_C;
        else
            Lbar_C = NaN;
        end

        volumeFromCycles = full(sum(sub.VolumeFromCycles));
        volumeGap = volumeFromCycles - V_C;
        volumeRelError = abs(volumeGap) / max(1, abs(V_C));

        batchSamples = full(sum(sub.BatchSamples));
        batchTransitions = full(sum(sub.BatchTransitions));

        batchT_C = full(sum(sub.Batch_T_C));

        if batchT_C > tol
            batchLbar_C = V_C / batchT_C;
        else
            batchLbar_C = NaN;
        end

        batchVolumeFromCycles = full(sum(sub.BatchVolumeFromCycles));
        batchVolumeGap = batchVolumeFromCycles - V_C;
        batchVolumeRelError = abs(batchVolumeGap) / max(1, abs(V_C));

        rows(end+1, :) = { ...
            b, ...
            cumulativeSamples, ...
            cumulativeTransitions, ...
            observedCycles, ...
            T_C, ...
            Lbar_C, ...
            volumeFromCycles, ...
            volumeGap, ...
            volumeRelError, ...
            batchSamples, ...
            batchTransitions, ...
            batchT_C, ...
            batchLbar_C, ...
            batchVolumeFromCycles, ...
            batchVolumeGap, ...
            batchVolumeRelError ...
        }; %#ok<AGROW>
    end

    if isempty(rows)
        Trace = table();
        return;
    end

    Trace = cell2table(rows, ...
        'VariableNames', { ...
            'Batch', ...
            'CumulativeSamples', ...
            'CumulativeTransitions', ...
            'ObservedCycles', ...
            'T_C', ...
            'Lbar_C', ...
            'VolumeFromCycles', ...
            'VolumeGap', ...
            'VolumeRelError', ...
            'BatchSamples', ...
            'BatchTransitions', ...
            'Batch_T_C', ...
            'Batch_Lbar_C', ...
            'BatchVolumeFromCycles', ...
            'BatchVolumeGap', ...
            'BatchVolumeRelError' ...
        });

    diagnostics.NumBatches = height(Trace);

    if height(Trace) >= 2
        diagnostics.FinalRelChangeT_C = abs(Trace.T_C(end) - Trace.T_C(end-1)) / ...
            max(1, abs(Trace.T_C(end)));

        diagnostics.FinalRelChangeLbar_C = abs(Trace.Lbar_C(end) - Trace.Lbar_C(end-1)) / ...
            max(1, abs(Trace.Lbar_C(end)));
    end

    validT = Trace.Batch_T_C(isfinite(Trace.Batch_T_C));
    validL = Trace.Batch_Lbar_C(isfinite(Trace.Batch_Lbar_C));

    if numel(validT) >= 2
        diagnostics.MCSE_T_C = std(validT) / sqrt(numel(validT));
    end

    if numel(validL) >= 2
        diagnostics.MCSE_Lbar_C = std(validL) / sqrt(numel(validL));
    end
end


function components = getPositiveSupportComponents(C, tol)
% getPositiveSupportComponents
% Return components of the positive support with positive edge-flow volume.

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

    Convergence = struct();
    Convergence.ComponentTrace = table();
    Convergence.Trace = table();
    Convergence.Diagnostics = struct( ...
        'Enabled', false, ...
        'NumComponentTraceRows', 0, ...
        'NumBatches', 0, ...
        'FinalRelChangeT_C', NaN, ...
        'FinalRelChangeLbar_C', NaN, ...
        'MCSE_T_C', NaN, ...
        'MCSE_Lbar_C', NaN, ...
        'StopReason', "empty" ...
    );

    MC = struct();
    MC.CycleTable = CycleTable;
    MC.TopCycles = table();
    MC.Estimates = Estimates;
    MC.ComponentTable = table();
    MC.Convergence = Convergence;
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
    opts.BatchTransitions = [];
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

            case 'batchtransitions'
                opts.BatchTransitions = value;

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

    if ~isempty(opts.BatchTransitions)
        validatePositiveIntegerLike(opts.BatchTransitions, 'BatchTransitions');

        if lower(string(opts.BudgetMode)) ~= "transitions"
            error('BatchTransitions currently requires BudgetMode=''transitions''.');
        end
    end

    if ~isinf(opts.MaxTransitions)
        validatePositiveIntegerLike(opts.MaxTransitions, 'MaxTransitions');
    end
end


function validatePositiveIntegerLike(x, name)

    if ~isnumeric(x) || ~isscalar(x) || ~isfinite(x) || x < 1 || x ~= round(x)
        error('%s must be a positive integer scalar.', name);
    end
end