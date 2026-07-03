function [T, status] = cdfdExactAllCycles(C, varargin)
% cdfdExactAllCycles
% Exact probabilistic cycle coefficients for all simple directed cycles.
%
% CONTRACT
%   This function returns exact all-cycle results only.
%   If full enumeration is infeasible or exceeds MaxCycles, the function
%   errors and returns no substantive result.
%
% INPUT
%   C   n x n non-negative circular component, typically C_BFF.
%
% OPTIONS
%   'Tol'                  default 1e-12
%   'MaxCycles'            default 50000
%   'CheckBalance'         default true
%   'CheckReconstruction'  default true
%   'MaxTotalCycleLength'  default 5e6
%   'MaxRuntimeSeconds'    default 60
%
% OUTPUT
%   T      table with columns:
%          Component, Cycle, Length, q, lambda, S
%
%   status struct with fields:
%          Complete, NumCycles, ReconstructionError, ComponentStatus
%          MaxCycles, MaxTotalCycleLength, MaxRuntimeSeconds,
%          TotalCycleLength, RuntimeSeconds
%
% NOTES
%   Requires MATLAB R2021a or later for allcycles when non-self cycles
%   of length at least 2 are present.
%   Self-loops are included as length-one simple cycles.
%   Non-self simple cycles are enumerated using allcycles.
%   This function is intended for small networks and validation examples.
%   For large networks, use streaming Monte Carlo.

    opts = parseExactAllOptions(varargin{:});
    tol = opts.Tol;

    runStartTic = tic;
    totalCycleLength = 0;

    validateSquareNonnegative(C, tol);

    % Numerical support convention.
    C(abs(C) <= tol) = 0;

    if opts.CheckBalance
        balanceResidual = full(max(abs(sum(C, 2) - sum(C, 1).')));
        scale = full(max(1, max(sum(C, 2))));

        if balanceResidual > 100 * tol * scale
            error(['Input C does not appear to be balanced. ' ...
                   'Max balance residual = %.3g. This function expects a circular component.'], ...
                   balanceResidual);
        end
    end

    n = size(C, 1);

    % Positive off-diagonal support for SCC decomposition and allcycles.
    % Self-loops are handled manually as length-one cycles.
    A = C > tol;
    A(1:n+1:end) = false;

    Gfull = digraph(A);
    comp = conncomp(Gfull, 'Type', 'strong');
    compIds = unique(comp);

    rowsComponent = zeros(0, 1);
    rowsCycle = cell(0, 1);
    rowsLength = zeros(0, 1);
    rowsQ = zeros(0, 1);
    rowsLambda = zeros(0, 1);
    rowsS = zeros(0, 1);

    compRows = {};

    status = struct();
    status.Complete = true;
    status.NumCycles = 0;
    status.ReconstructionError = NaN;
    status.ComponentStatus = table();
    status.MaxCycles = opts.MaxCycles;
    status.MaxTotalCycleLength = opts.MaxTotalCycleLength;
    status.MaxRuntimeSeconds = opts.MaxRuntimeSeconds;
    status.TotalCycleLength = 0;
    status.RuntimeSeconds = NaN;

    hasAllCycles = ~isempty(which('allcycles'));

    for idx = 1:numel(compIds)
        cc = compIds(idx);
        nodes = find(comp == cc);

        Csub = C(nodes, nodes);
        m = numel(nodes);

        Asub = Csub > tol;
        Asub(1:m+1:end) = false;

        selfLocal = find(diag(Csub) > tol).';

        if nnz(Asub) == 0 && isempty(selfLocal)
            continue;
        end

        cyclesLocal = cell(0, 1);

        % Add self-loops manually as length-one cycles.
        for a = selfLocal
            cyclesLocal{end+1, 1} = a; %#ok<AGROW>
        end

        remaining = opts.MaxCycles - status.NumCycles;

        if numel(cyclesLocal) > remaining
            error(['Full cycle enumeration exceeded MaxCycles=%d. ' ...
                   'No exact all-cycle result returned. Use streaming Monte Carlo ' ...
                   'or increase MaxCycles if full enumeration is intended.'], ...
                   opts.MaxCycles);
        end

        remainingAfterSelf = remaining - numel(cyclesLocal);

        componentCycleLength = sum(cellfun(@numel, cyclesLocal));

        if totalCycleLength + componentCycleLength > opts.MaxTotalCycleLength
            error(['Exact cycle enumeration exceeded MaxTotalCycleLength=%d. ' ...
                'No exact all-cycle result returned. Use streaming Monte Carlo ' ...
                'or increase MaxTotalCycleLength if full enumeration is intended.'], ...
                opts.MaxTotalCycleLength);
        end

        % Enumerate non-self simple cycles using allcycles.
        if nnz(Asub) > 0
            if ~hasAllCycles
                error(['MATLAB allcycles was not found. It was introduced in R2021a. ' ...
                       'Use MATLAB R2021a or later, or provide an alternative simple-cycle enumerator.']);
            end

            if remainingAfterSelf <= 0
                error(['Full cycle enumeration exceeded MaxCycles=%d. ' ...
                       'No exact all-cycle result returned. Use streaming Monte Carlo ' ...
                       'or increase MaxCycles if full enumeration is intended.'], ...
                       opts.MaxCycles);
            end

            G = digraph(Asub);

            % Request one extra cycle to detect whether this component exceeds
            % the remaining all-cycle budget.
            cyclesNonself = guardedAllcycles(G, remainingAfterSelf, opts);

            nonselfCycleLength = sum(cellfun(@numel, cyclesNonself));

            if totalCycleLength + componentCycleLength + nonselfCycleLength > opts.MaxTotalCycleLength
                error(['Exact cycle enumeration exceeded MaxTotalCycleLength=%d. ' ...
                    'No exact all-cycle result returned. Use streaming Monte Carlo ' ...
                    'or increase MaxTotalCycleLength if full enumeration is intended.'], ...
                    opts.MaxTotalCycleLength);
            end

            cyclesLocal = [cyclesLocal; cyclesNonself(:)];

            %componentCycleLength = sum(cellfun(@numel, cyclesLocal));
            %totalCycleLength = totalCycleLength + componentCycleLength;

        end

        componentCycleLength = sum(cellfun(@numel, cyclesLocal));

        if totalCycleLength + componentCycleLength > opts.MaxTotalCycleLength
            error(['Exact cycle enumeration exceeded MaxTotalCycleLength=%d. ' ...
                'No exact all-cycle result returned. Use streaming Monte Carlo ' ...
                'or increase MaxTotalCycleLength if full enumeration is intended.'], ...
                opts.MaxTotalCycleLength);
        end

        totalCycleLength = totalCycleLength + componentCycleLength;

        for r = 1:numel(cyclesLocal)

            guardExactRuntime(runStartTic, opts);
            cycLocal = cyclesLocal{r}(:).';
            cycGlobal = nodes(cycLocal);

            [q, lambda, ~, S] = qqWeightLocal(Csub, cycLocal, tol);

            rowsComponent(end+1, 1) = cc; %#ok<AGROW>
            rowsCycle{end+1, 1} = cycGlobal; %#ok<AGROW>
            rowsLength(end+1, 1) = numel(cycGlobal); %#ok<AGROW>
            rowsQ(end+1, 1) = q; %#ok<AGROW>
            rowsLambda(end+1, 1) = lambda; %#ok<AGROW>
            rowsS(end+1, 1) = S; %#ok<AGROW>

            status.NumCycles = status.NumCycles + 1;
        end

        compRows(end+1, :) = { ...
            cc, ...
            m, ...
            nnz(Csub > tol), ...
            nnz(Asub), ...
            numel(selfLocal), ...
            numel(cyclesLocal) ...
        }; %#ok<AGROW>
    end

    T = table( ...
        rowsComponent, ...
        rowsCycle, ...
        rowsLength, ...
        rowsQ, ...
        rowsLambda, ...
        rowsS, ...
        'VariableNames', {'Component', 'Cycle', 'Length', 'q', 'lambda', 'S'} ...
    );

    if ~isempty(compRows)
        status.ComponentStatus = cell2table(compRows, ...
            'VariableNames', { ...
                'Component', ...
                'NumNodes', ...
                'NumEdges', ...
                'NumNonSelfEdges', ...
                'NumSelfLoops', ...
                'NumCyclesFound' ...
            });
    end

    if opts.CheckReconstruction
        status.ReconstructionError = reconstructionErrorFromTable(C, T);

        scale = full(max(1, norm(C, 'fro')));
        if status.ReconstructionError > 1e-8 * scale
            error(['Exact cycle reconstruction failed. Frobenius error = %.3g. ' ...
                   'Check that C is a balanced circular component, that all supported cycles ' ...
                   'were enumerated, and that Tol is appropriate.'], ...
                   status.ReconstructionError);
        end
    end

    status.TotalCycleLength = totalCycleLength;
    status.RuntimeSeconds = toc(runStartTic);

end


function [q, lambda, logq, S] = qqWeightLocal(Csub, cycleLocal, tol)
% qqWeightLocal
% Qian--Qian exact coefficient for one cycle in one SCC.
%
% cycleLocal is a row vector of local node indices. Self-loops are
% represented by a length-one vector, e.g. [i].

    Csub = full(Csub);

    m = size(Csub, 1);
    s = full(sum(Csub, 2));
    S = full(sum(s));

    if S <= tol
        error('Component has zero circular edge-flow volume.');
    end

    if any(s <= tol)
        error('Component contains a node with zero outflow.');
    end

    P = Csub ./ s;
    L = eye(m) - P;

    closed = [cycleLocal cycleLocal(1)];

    logProdP = 0;
    for k = 1:numel(cycleLocal)
        i = closed(k);
        j = closed(k+1);
        pij = P(i, j);

        if pij <= tol
            error('Cycle transition has zero probability inside component.');
        end

        logProdP = logProdP + log(pij);
    end

    cycleSet = unique(cycleLocal);
    complement = setdiff(1:m, cycleSet);

    logMinorNumerator = logdetAbs(L(complement, complement));

    logDenTerms = -inf(1, m);
    for k = 1:m
        idx = setdiff(1:m, k);
        logDenTerms(k) = logdetAbs(L(idx, idx));
    end

    logDenominator = logsumexp(logDenTerms);

    logq = logProdP + logMinorNumerator - logDenominator;
    q = exp(logq);
    lambda = S * q;
end


function err = reconstructionErrorFromTable(C, T)
% reconstructionErrorFromTable
% Reconstruct C from cycle table and return Frobenius error.

    n = size(C, 1);
    Chat = zeros(n);

    for r = 1:height(T)
        cyc = T.Cycle{r};
        lam = full(T.lambda(r));

        closed = [cyc(:).' cyc(1)];

        for k = 1:numel(cyc)
            i = closed(k);
            j = closed(k+1);
            Chat(i, j) = Chat(i, j) + lam;
        end
    end

    err = norm(Chat - C, 'fro');
end


function ld = logdetAbs(A)
% logdetAbs
% Log absolute determinant, with det(empty)=1.

    if isempty(A)
        ld = 0;
        return;
    end

    [~, U, ~] = lu(A, 'vector');
    diagU = diag(U);

    if any(abs(diagU) == 0)
        ld = -inf;
    else
        ld = sum(log(abs(diagU)));
    end
end


function y = logsumexp(x)
% logsumexp
% Stable log(sum(exp(x))).

    xmax = max(x);

    if isinf(xmax)
        y = xmax;
    else
        y = xmax + log(sum(exp(x - xmax)));
    end
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


function opts = parseExactAllOptions(varargin)

    opts = struct();
    opts.Tol = 1e-12;
    opts.MaxCycles = 50000;
    opts.MaxTotalCycleLength = 5e6;
    opts.MaxRuntimeSeconds = 60;
    opts.CheckBalance = true;
    opts.CheckReconstruction = true;

    if mod(numel(varargin), 2) ~= 0
        error('Options must be name-value pairs.');
    end

    for k = 1:2:numel(varargin)
        name = varargin{k};
        value = varargin{k+1};

        switch lower(name)
            case 'tol'
                opts.Tol = value;
            case 'maxcycles'
                opts.MaxCycles = value;
            case 'checkbalance'
                opts.CheckBalance = value;
            case 'checkreconstruction'
                opts.CheckReconstruction = value;
            case 'maxtotalcyclelength'
                opts.MaxTotalCycleLength = value;
            case 'maxruntimeseconds'
                opts.MaxRuntimeSeconds = value;
            otherwise
                error('Unknown option: %s', name);
        end
    end

    validatePositiveIntegerLike(opts.MaxCycles, 'MaxCycles');
    validatePositiveIntegerLike(opts.MaxTotalCycleLength, 'MaxTotalCycleLength');
    validatePositiveFiniteScalar(opts.MaxRuntimeSeconds, 'MaxRuntimeSeconds');
    validatePositiveFiniteScalar(opts.Tol, 'Tol');

    if ~islogical(opts.CheckBalance) && ~ismember(opts.CheckBalance, [0 1])
        error('CheckBalance must be true or false.');
    end

    if ~islogical(opts.CheckReconstruction) && ~ismember(opts.CheckReconstruction, [0 1])
        error('CheckReconstruction must be true or false.');
    end

end

function cycles = guardedAllcycles(G, remainingBudget, opts)
% guardedAllcycles
% Hard-guarded wrapper around MATLAB allcycles.
%
% This intentionally does not fall back to unguarded allcycles. Exact
% all-cycle enumeration must fail rather than launch an unbounded
% enumeration.

    validatePositiveIntegerLike(remainingBudget, 'remainingBudget');

    if remainingBudget <= 0
        error(['Full cycle enumeration exceeded MaxCycles=%d. ' ...
               'No exact all-cycle result returned. Use streaming Monte Carlo ' ...
               'or increase MaxCycles if full enumeration is intended.'], ...
               opts.MaxCycles);
    end

    cycles = allcycles(G, ...
        'MinCycleLength', 2, ...
        'MaxNumCycles', remainingBudget + 1);

    if numel(cycles) > remainingBudget
        error(['Full cycle enumeration exceeded MaxCycles=%d. ' ...
               'No exact all-cycle result returned. Use streaming Monte Carlo ' ...
               'or increase MaxCycles if full enumeration is intended.'], ...
               opts.MaxCycles);
    end
end


function guardExactRuntime(runStartTic, opts)
% guardExactRuntime
% Runtime guard for post-enumeration work.
%
% Note: this cannot interrupt MATLAB while it is inside allcycles. The
% hard guard against excessive cycle materialisation is MaxNumCycles.

    elapsed = toc(runStartTic);

    if elapsed > opts.MaxRuntimeSeconds
        error(['Exact cycle analysis exceeded MaxRuntimeSeconds=%.3g. ' ...
               'No exact all-cycle result returned. Use streaming Monte Carlo ' ...
               'or increase MaxRuntimeSeconds if full exact weighting is intended.'], ...
               opts.MaxRuntimeSeconds);
    end
end


function validatePositiveIntegerLike(x, name)
    if ~isnumeric(x) || ~isscalar(x) || ~isfinite(x) || x < 0 || floor(x) ~= x
        error('%s must be a non-negative finite integer scalar.', name);
    end
end


function validatePositiveFiniteScalar(x, name)
    if ~isnumeric(x) || ~isscalar(x) || ~isfinite(x) || x <= 0
        error('%s must be a positive finite scalar.', name);
    end
end

