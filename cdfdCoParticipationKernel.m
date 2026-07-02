function Kout = cdfdCoParticipationKernel(Tcycles, n, varargin)
% cdfdCoParticipationKernel
% Edge co-participation kernel from probabilistic cycle coefficients.
%
% PURPOSE
%   Given a cycle table with coefficients lambda_gamma, compute the
%   edge-by-edge co-participation kernel
%
%       K_{ef} = sum_gamma lambda_gamma * 1{e in gamma} * 1{f in gamma}.
%
%   The diagonal satisfies
%
%       K_{ee} = c_e,
%
%   where c_e is the circular flow on edge e reconstructed from the cycle
%   table. The normalised kernel is
%
%       Khat_{ef} = K_{ef} / sqrt(c_e c_f),
%
%   with zero entries when c_e or c_f is numerically zero.
%
% INPUT
%   Tcycles   table with at least columns:
%               Cycle, lambda
%
%             Usually returned by cdfdExactAllCycles(C).
%
%   n         optional number of nodes. If omitted or empty, n is inferred
%             from the maximum node index appearing in Tcycles.
%
% OPTIONS
%   'Tol'              default 1e-12
%   'WeightVariable'   default 'lambda'
%   'CycleVariable'    default 'Cycle'
%   'MakePairTable'    default false. If true, create a long-form table of
%                      edge-pair co-participation values.
%   'IncludeZeroPairs' default false. Only used when MakePairTable is true.
%
% OUTPUT
%   Kout      struct with fields:
%
%             K                  raw edge co-participation matrix
%             Khat               normalised co-participation matrix
%             B                  sparse edge-cycle incidence matrix
%             EdgeTable          table mapping edge index to From, To, Cedge
%             EdgePairTable      optional long-form edge-pair table
%             Cedge              vector of reconstructed edge circular flows
%             CedgeMatrix        n x n reconstructed circular component
%             Lambda             cycle coefficient vector
%             NumNodes           number of nodes
%             NumEdges           number of edges appearing in listed cycles
%             NumCycles          number of cycles
%             Options            parsed options
%
% NOTES
%   - Self-loops are allowed and are treated as length-one cycles.
%   - A cycle [i] contributes to edge i -> i.
%   - Cycles may be supplied open, e.g. [1 2 3], or closed, e.g. [1 2 3 1].
%   - If Tcycles is a partial/estimated cycle table, K is the
%     co-participation kernel for the listed/estimated cycles only.

    if nargin < 2
        n = [];
        varargin = {};
    elseif ischar(n) || isstring(n)
        varargin = [{n}, varargin];
        n = [];
    end

    opts = parseKernelOptions(varargin{:});
    tol = opts.Tol;

    validateCycleTable(Tcycles, opts);

    cycles = Tcycles.(opts.CycleVariable);
    lambda = full(Tcycles.(opts.WeightVariable));
    lambda = lambda(:);

    if any(~isfinite(lambda)) || any(lambda < -tol)
        error('Cycle weights must be finite and non-negative up to tolerance.');
    end

    lambda(abs(lambda) <= tol) = 0;

    numCycles = height(Tcycles);

    if isempty(n)
        n = inferNumNodes(cycles);
    end

    validateNumNodes(n);

    % First pass: normalise cycles and collect directed edges.
    normCycles = cell(numCycles, 1);
    edgeRows = zeros(0, 2);

    for r = 1:numCycles
        cyc = normaliseCycle(cycles{r});

        if isempty(cyc)
            error('Cycle table contains an empty cycle.');
        end

        validateCycleNodes(cyc, n, r);

        normCycles{r} = cyc;

        closed = [cyc cyc(1)];

        for k = 1:numel(cyc)
            i = closed(k);
            j = closed(k+1);
            edgeRows(end+1, :) = [i j]; %#ok<AGROW>
        end
    end

    if isempty(edgeRows)
        Kout = emptyOutput(n, numCycles, lambda, opts);
        return;
    end

    % Unique directed edge list.
    edgeRows = unique(edgeRows, 'rows', 'stable');
    numEdges = size(edgeRows, 1);

    edgeKey = containers.Map('KeyType', 'char', 'ValueType', 'double');

    for e = 1:numEdges
        edgeKey(edgeToKey(edgeRows(e, 1), edgeRows(e, 2))) = e;
    end

    % Sparse edge-cycle incidence matrix B.
    % B(e,gamma)=1 iff edge e appears in cycle gamma.
    ii = zeros(0, 1);
    jj = zeros(0, 1);
    vv = zeros(0, 1);

    for r = 1:numCycles
        cyc = normCycles{r};
        closed = [cyc cyc(1)];

        edgeIdxInCycle = zeros(numel(cyc), 1);

        for k = 1:numel(cyc)
            i = closed(k);
            j = closed(k+1);

            key = edgeToKey(i, j);
            edgeIdxInCycle(k) = edgeKey(key);
        end

        % A simple cycle should not repeat directed edges, but unique keeps
        % the incidence interpretation safe.
        edgeIdxInCycle = unique(edgeIdxInCycle);

        ii = [ii; edgeIdxInCycle(:)]; %#ok<AGROW>
        jj = [jj; r * ones(numel(edgeIdxInCycle), 1)]; %#ok<AGROW>
        vv = [vv; ones(numel(edgeIdxInCycle), 1)]; %#ok<AGROW>
    end

    B = sparse(ii, jj, vv, numEdges, numCycles);

    % K = B diag(lambda) B'.
    K = B * spdiags(lambda, 0, numCycles, numCycles) * B.';
    K = full(K);
    K(abs(K) <= tol) = 0;

    Cedge = full(diag(K));

    denom = sqrt(Cedge * Cedge.');
    Khat = zeros(size(K));

    mask = denom > tol;
    Khat(mask) = K(mask) ./ denom(mask);
    Khat(abs(Khat) <= tol) = 0;

    CedgeMatrix = zeros(n);

    for e = 1:numEdges
        i = edgeRows(e, 1);
        j = edgeRows(e, 2);
        CedgeMatrix(i, j) = Cedge(e);
    end

    From = edgeRows(:, 1);
    To = edgeRows(:, 2);
    EdgeIndex = (1:numEdges).';
    IsSelfLoop = From == To;

    EdgeTable = table( ...
        EdgeIndex, ...
        From, ...
        To, ...
        Cedge, ...
        IsSelfLoop, ...
        'VariableNames', { ...
            'EdgeIndex', ...
            'From', ...
            'To', ...
            'Cedge', ...
            'IsSelfLoop' ...
        } ...
    );

    Kout = struct();
    Kout.K = K;
    Kout.Khat = Khat;
    Kout.B = B;
    Kout.EdgeTable = EdgeTable;
    Kout.Cedge = Cedge;
    Kout.CedgeMatrix = CedgeMatrix;
    Kout.Lambda = lambda;
    Kout.NumNodes = n;
    Kout.NumEdges = numEdges;
    Kout.NumCycles = numCycles;
    Kout.Options = opts;

    if opts.MakePairTable
        Kout.EdgePairTable = makeEdgePairTable(K, Khat, EdgeTable, opts.IncludeZeroPairs, tol);
    else
        Kout.EdgePairTable = table();
    end
end


function pairTable = makeEdgePairTable(K, Khat, EdgeTable, includeZeroPairs, tol)
% makeEdgePairTable
% Long-form edge-pair table for co-participation values.

    m = height(EdgeTable);

    EdgeIndex1 = zeros(0, 1);
    EdgeIndex2 = zeros(0, 1);
    From1 = zeros(0, 1);
    To1 = zeros(0, 1);
    From2 = zeros(0, 1);
    To2 = zeros(0, 1);
    CoParticipation = zeros(0, 1);
    NormalisedCoParticipation = zeros(0, 1);

    for a = 1:m
        for b = 1:m
            if ~includeZeroPairs && abs(K(a, b)) <= tol
                continue;
            end

            EdgeIndex1(end+1, 1) = EdgeTable.EdgeIndex(a); %#ok<AGROW>
            EdgeIndex2(end+1, 1) = EdgeTable.EdgeIndex(b); %#ok<AGROW>

            From1(end+1, 1) = EdgeTable.From(a); %#ok<AGROW>
            To1(end+1, 1) = EdgeTable.To(a); %#ok<AGROW>
            From2(end+1, 1) = EdgeTable.From(b); %#ok<AGROW>
            To2(end+1, 1) = EdgeTable.To(b); %#ok<AGROW>

            CoParticipation(end+1, 1) = K(a, b); %#ok<AGROW>
            NormalisedCoParticipation(end+1, 1) = Khat(a, b); %#ok<AGROW>
        end
    end

    pairTable = table( ...
        EdgeIndex1, ...
        EdgeIndex2, ...
        From1, ...
        To1, ...
        From2, ...
        To2, ...
        CoParticipation, ...
        NormalisedCoParticipation ...
    );
end


function Kout = emptyOutput(n, numCycles, lambda, opts)
% emptyOutput
% Return a consistent empty output struct.

    Kout = struct();
    Kout.K = zeros(0);
    Kout.Khat = zeros(0);
    Kout.B = sparse(0, numCycles);
    Kout.EdgeTable = table();
    Kout.EdgePairTable = table();
    Kout.Cedge = zeros(0, 1);
    Kout.CedgeMatrix = zeros(n);
    Kout.Lambda = lambda;
    Kout.NumNodes = n;
    Kout.NumEdges = 0;
    Kout.NumCycles = numCycles;
    Kout.Options = opts;
end


function cyc = normaliseCycle(cyc)
% normaliseCycle
% Convert cycle to open representation.
%
% Examples:
%   [1 2 3 1] -> [1 2 3]
%   [1 1]     -> [1]
%   [1]       -> [1]

    cyc = cyc(:).';

    if numel(cyc) >= 2 && cyc(end) == cyc(1)
        cyc = cyc(1:end-1);
    end
end


function key = edgeToKey(i, j)

    key = sprintf('%d_%d', i, j);
end


function n = inferNumNodes(cycles)

    n = 0;

    for r = 1:numel(cycles)
        cyc = normaliseCycle(cycles{r});

        if ~isempty(cyc)
            n = max(n, max(cyc));
        end
    end
end


function validateCycleNodes(cyc, n, row)

    if any(cyc < 1) || any(cyc > n) || any(cyc ~= round(cyc))
        error('Cycle table contains invalid node index in row %d.', row);
    end

    if numel(unique(cyc)) ~= numel(cyc)
        error('Cycle table contains a non-simple cycle or repeated node in row %d.', row);
    end
end


function validateNumNodes(n)

    if ~isnumeric(n) || ~isscalar(n) || ~isfinite(n) || n < 0 || n ~= round(n)
        error('n must be a non-negative integer scalar.');
    end
end


function validateCycleTable(Tcycles, opts)

    if ~istable(Tcycles)
        error('Tcycles must be a MATLAB table.');
    end

    required = {opts.CycleVariable, opts.WeightVariable};

    for k = 1:numel(required)
        if ~ismember(required{k}, Tcycles.Properties.VariableNames)
            error('Tcycles is missing required variable: %s.', required{k});
        end
    end

    cycles = Tcycles.(opts.CycleVariable);

    if ~iscell(cycles)
        error('Tcycles.%s must be a cell array of node-index vectors.', opts.CycleVariable);
    end
end


function opts = parseKernelOptions(varargin)

    opts = struct();
    opts.Tol = 1e-12;
    opts.WeightVariable = 'lambda';
    opts.CycleVariable = 'Cycle';
    opts.MakePairTable = false;
    opts.IncludeZeroPairs = false;

    if mod(numel(varargin), 2) ~= 0
        error('Options must be name-value pairs.');
    end

    for k = 1:2:numel(varargin)
        name = varargin{k};
        value = varargin{k+1};

        switch lower(name)
            case 'tol'
                opts.Tol = value;

            case 'weightvariable'
                opts.WeightVariable = value;

            case 'cyclevariable'
                opts.CycleVariable = value;

            case 'makepairtable'
                opts.MakePairTable = value;

            case 'includezeropairs'
                opts.IncludeZeroPairs = value;

            otherwise
                error('Unknown option: %s.', name);
        end
    end
end