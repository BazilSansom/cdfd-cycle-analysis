function Top = cdfdTopCycles(Tcycles, varargin)
% cdfdTopCycles
% Ranked table of important cycles from probabilistic cycle coefficients.
%
% PURPOSE
%   Converts a cycle table, usually returned by cdfdExactAllCycles(C), into
%   a ranked table of feedback-loop importance measures:
%
%       I_thr(gamma) = lambda_gamma
%       I_vol(gamma) = lambda_gamma * |gamma|
%
%   The output includes throughput shares, volume shares, and cumulative
%   coverage after ranking.
%
% INPUT
%   Tcycles   table with at least columns:
%               Cycle, Length, lambda
%
% OPTIONS
%   'TopK'             default Inf. Number of ranked cycles to return.
%   'RankBy'           default 'throughput'. One of:
%                        'throughput'  rank by lambda
%                        'volume'      rank by lambda * Length
%                        'length'      rank by Length
%   'Descending'       default true.
%   'Tol'              default 1e-12.
%   'WeightVariable'   default 'lambda'.
%   'LengthVariable'   default 'Length'.
%   'CycleVariable'    default 'Cycle'.
%
% OUTPUT
%   Top   table with columns:
%           Rank
%           CycleIndex
%           Component             if present in Tcycles
%           Cycle
%           Label
%           Length
%           lambda
%           VolumeContribution
%           ThroughputShare
%           VolumeShare
%           CumulativeThroughputShare
%           CumulativeVolumeShare
%           RankBy
%           RankValue
%           Polarity              if present in Tcycles
%           IsPositive            if present in Tcycles
%           IsNegative            if present in Tcycles
%
% NOTES
%   - Self-loops are allowed and labelled as "i -> i".
%   - Two-cycles are labelled as "i <-> j".
%   - Longer cycles are labelled with closure, e.g. "1 -> 2 -> 3 -> 1".
%   - If Tcycles is a partial/estimated table, coverage shares are relative
%     to the listed/estimated cycles only.

    opts = parseTopCycleOptions(varargin{:});
    tol = opts.Tol;

    validateCycleTable(Tcycles, opts);

    numCycles = height(Tcycles);

    if numCycles == 0
        Top = table();
        return;
    end

    cycles = Tcycles.(opts.CycleVariable);
    lambda = full(Tcycles.(opts.WeightVariable));
    lengths = full(Tcycles.(opts.LengthVariable));

    lambda = lambda(:);
    lengths = lengths(:);

    if any(~isfinite(lambda)) || any(lambda < -tol)
        error('Cycle weights must be finite and non-negative up to tolerance.');
    end

    if any(~isfinite(lengths)) || any(lengths < 1) || any(abs(lengths - round(lengths)) > tol)
        error('Cycle lengths must be positive integers.');
    end

    lambda(abs(lambda) <= tol) = 0;
    lengths = round(lengths);

    labels = strings(numCycles, 1);
    normalisedCycles = cell(numCycles, 1);

    for r = 1:numCycles
        cyc = normaliseCycle(cycles{r});

        if isempty(cyc)
            error('Cycle table contains an empty cycle in row %d.', r);
        end

        if numel(cyc) ~= lengths(r)
            error('Cycle length mismatch in row %d.', r);
        end

        if numel(unique(cyc)) ~= numel(cyc)
            error('Cycle table contains a non-simple cycle or repeated node in row %d.', r);
        end

        normalisedCycles{r} = cyc;
        labels(r) = makeCycleLabel(cyc);
    end

    volumeContribution = lambda .* lengths;

    T_C = full(sum(lambda));
    V_C = full(sum(volumeContribution));

    throughputShare = safeShareVector(lambda, T_C, tol);
    volumeShare = safeShareVector(volumeContribution, V_C, tol);

    rankBy = lower(string(opts.RankBy));

    switch rankBy
        case "throughput"
            rankValue = lambda;

        case "volume"
            rankValue = volumeContribution;

        case "length"
            rankValue = lengths;

        otherwise
            error('Unknown RankBy option: %s.', opts.RankBy);
    end

    if opts.Descending
        [~, order] = sort(rankValue, 'descend');
    else
        [~, order] = sort(rankValue, 'ascend');
    end

    topK = opts.TopK;

    if isempty(topK) || isinf(topK)
        topK = numCycles;
    end

    if ~isnumeric(topK) || ~isscalar(topK) || topK < 0 || ~isfinite(topK)
        error('TopK must be a non-negative scalar or Inf.');
    end

    topK = min(numCycles, floor(topK));
    order = order(1:topK);

    Rank = (1:topK).';
    CycleIndex = order(:);

    Cycle = normalisedCycles(order);
    Label = labels(order);
    Length = lengths(order);
    lambdaOut = lambda(order);
    VolumeContribution = volumeContribution(order);
    ThroughputShare = throughputShare(order);
    VolumeShare = volumeShare(order);

    CumulativeThroughputShare = cumsum(ThroughputShare);
    CumulativeVolumeShare = cumsum(VolumeShare);

    RankBy = repmat(rankBy, topK, 1);
    RankValue = rankValue(order);

    Top = table( ...
        Rank, ...
        CycleIndex, ...
        Cycle, ...
        Label, ...
        Length, ...
        lambdaOut, ...
        VolumeContribution, ...
        ThroughputShare, ...
        VolumeShare, ...
        CumulativeThroughputShare, ...
        CumulativeVolumeShare, ...
        RankBy, ...
        RankValue, ...
        'VariableNames', { ...
            'Rank', ...
            'CycleIndex', ...
            'Cycle', ...
            'Label', ...
            'Length', ...
            opts.WeightVariable, ...
            'VolumeContribution', ...
            'ThroughputShare', ...
            'VolumeShare', ...
            'CumulativeThroughputShare', ...
            'CumulativeVolumeShare', ...
            'RankBy', ...
            'RankValue' ...
        } ...
    );

    % Add optional variables from the original cycle table when available.
    if ismember('Component', Tcycles.Properties.VariableNames)
        Top = addvars(Top, Tcycles.Component(order), ...
            'Before', 'Cycle', ...
            'NewVariableNames', 'Component');
    end

    optionalVars = {'Polarity', 'IsPositive', 'IsNegative'};

    for k = 1:numel(optionalVars)
        v = optionalVars{k};

        if ismember(v, Tcycles.Properties.VariableNames)
            Top.(v) = Tcycles.(v)(order);
        end
    end
end


function s = safeShareVector(x, total, tol)

    if total > tol
        s = x ./ total;
    else
        s = NaN(size(x));
    end
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


function label = makeCycleLabel(cyc)

    if numel(cyc) == 1
        label = sprintf('%d -> %d', cyc(1), cyc(1));

    elseif numel(cyc) == 2
        label = sprintf('%d <-> %d', cyc(1), cyc(2));

    else
        closed = [cyc(:).' cyc(1)];
        label = strjoin(string(closed), " -> ");
    end

    label = string(label);
end


function validateCycleTable(Tcycles, opts)

    if ~istable(Tcycles)
        error('Tcycles must be a MATLAB table.');
    end

    required = {opts.CycleVariable, opts.LengthVariable, opts.WeightVariable};

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


function opts = parseTopCycleOptions(varargin)

    opts = struct();
    opts.TopK = Inf;
    opts.RankBy = 'throughput';
    opts.Descending = true;
    opts.Tol = 1e-12;
    opts.WeightVariable = 'lambda';
    opts.LengthVariable = 'Length';
    opts.CycleVariable = 'Cycle';

    if mod(numel(varargin), 2) ~= 0
        error('Options must be name-value pairs.');
    end

    for k = 1:2:numel(varargin)
        name = varargin{k};
        value = varargin{k+1};

        switch lower(name)
            case 'topk'
                opts.TopK = value;

            case 'rankby'
                opts.RankBy = value;

            case 'descending'
                opts.Descending = value;

            case 'tol'
                opts.Tol = value;

            case 'weightvariable'
                opts.WeightVariable = value;

            case 'lengthvariable'
                opts.LengthVariable = value;

            case 'cyclevariable'
                opts.CycleVariable = value;

            otherwise
                error('Unknown option: %s.', name);
        end
    end
end