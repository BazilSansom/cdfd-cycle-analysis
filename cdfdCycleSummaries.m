function [S, diagnostics] = cdfdCycleSummaries(W, C, Tcycles, varargin)
% cdfdCycleSummaries
% Network-level summaries from a CDFD circular component and cycle table.
%
% PURPOSE
%   Aggregates exact or estimated probabilistic cycle coefficients into the
%   paper-level quantities:
%
%       V_W, V_C, T_C, Lbar_C, D_net, Gamma_vol, Gamma_thr, R_thr,
%       N_eff, length summaries, and top-k coverage tables.
%
% INPUT
%   W         n x n non-negative original weighted directed network.
%   C         n x n non-negative circular component, typically C_BFF.
%   Tcycles   table with at least columns:
%             Cycle, Length, lambda
%
%             Usually returned by cdfdExactAllCycles(C).
%
% OPTIONS
%   'Tol'                     default 1e-12
%   'CheckVolumeIdentity'      default true
%   'CheckReconstruction'      default true
%   'ErrorOnCheckFailure'      default false
%   'TopK'                    default [1 5 10 25 50 100]
%   'WeightVariable'           default 'lambda'
%   'LengthVariable'           default 'Length'
%   'CycleVariable'            default 'Cycle'
%
% OUTPUT
%   S             struct of network-level summary quantities, with fields:
%
%                 Basic counts:
%                   NumNodes                  number of nodes
%                   NumEdges                  number of positive edges in W
%                   NumCircularEdges          number of positive edges in C
%                   NumCycles                 number of cycles in Tcycles
%                   NumSelfLoops              number of length-one cycles
%                   NumComponentsInCycleTable number of SCCs represented in Tcycles
%
%                 Edge-flow volumes:
%                   V_W                       total edge-flow volume, sum(W(:))
%                   V_C                       circular edge-flow volume, sum(C(:))
%                   VolumeFromCycles          sum(lambda .* Length)
%                   VolumeGap                 VolumeFromCycles - V_C
%                   VolumeAbsError            abs(VolumeGap)
%                   VolumeRelError            relative volume identity error
%
%                 Circular throughput:
%                   T_C                       completed circular throughput, sum(lambda)
%                   Lbar_C                    throughput-weighted mean cycle length,
%                                             V_C / T_C when T_C > 0
%
%                 Directional throughput:
%                   ImbalanceVector           b_i = out_i(W) - in_i(W)
%                   ImbalanceResidual         abs(sum(ImbalanceVector))
%                   D_net                     sum(max(b_i,0))
%                   D_net_abs_form            0.5 * sum(abs(b_i))
%                   D_net_gap                 D_net - D_net_abs_form
%                   D_net_abs_error           abs(D_net_gap)
%
%                 Circularity measures:
%                   Gamma_vol                 V_C / V_W when V_W > 0
%                   Gamma_thr                 T_C / (T_C + D_net), with zero-flow
%                                             conventions
%                   R_thr                     T_C / D_net when defined
%
%                 Cycle-length summaries:
%                   MinLength                 minimum cycle length
%                   MaxLength                 maximum cycle length
%                   UnweightedMeanLength      unweighted mean of listed cycle lengths
%                   ThroughputWeightedMeanLength
%                                             same as Lbar_C
%
%                 Cycle concentration:
%                   HHI                       Herfindahl concentration of throughput
%                                             shares lambda / T_C
%                   N_eff                     inverse-Herfindahl effective number
%                                             of throughput-carrying cycles
%
%                 Coverage tables:
%                   TopKByThroughput          top-k coverage table after ranking
%                                             cycles by lambda
%                   TopKByVolume              top-k coverage table after ranking
%                                             cycles by lambda .* Length
%
%   diagnostics   struct of consistency checks, with fields:
%                   BalanceResidualC
%                   BalanceResidualW
%                   ImbalanceSumResidual
%                   VolumeIdentityChecked
%                   VolumeIdentityAbsError
%                   VolumeIdentityRelError
%                   VolumeIdentityPass
%                   ReconstructionChecked
%                   ReconstructionError
%                   ReconstructionRelError
%                   ReconstructionPass
%                   DirectionalThroughputIdentityError
%                   DirectionalThroughputIdentityPass
%
% NOTES
%   This function is intentionally not a full workflow/orchestrator.
%   It does not compute BFF and does not choose between exact enumeration
%   and Monte Carlo. A future wrapper such as cdfdCycleAnalysis(W,...) can
%   call bff_decompose, choose the cycle method, and then call this utility.

    opts = parseCycleSummaryOptions(varargin{:});
    tol = opts.Tol;

    validateSquareFinite(W, 'W');
    validateSquareFinite(C, 'C');

    if size(W, 1) ~= size(C, 1)
        error('W and C must have the same dimensions.');
    end

    if any(W(:) < -tol)
        error('W must be non-negative up to numerical tolerance.');
    end

    if any(C(:) < -tol)
        error('C must be non-negative up to numerical tolerance.');
    end

    Wwork = W;
    Cwork = C;

    Wwork(abs(Wwork) <= tol) = 0;
    Cwork(abs(Cwork) <= tol) = 0;

    validateCycleTable(Tcycles, opts);
    
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

    n = size(Wwork, 1);

    S = struct();

    % Basic network sizes.
    S.NumNodes = n;
    S.NumEdges = nnz(Wwork > tol);
    S.NumCircularEdges = nnz(Cwork > tol);
    S.NumCycles = height(Tcycles);
    S.NumSelfLoops = sum(lengths == 1);

    if ismember('Component', Tcycles.Properties.VariableNames) && ~isempty(Tcycles)
        S.NumComponentsInCycleTable = numel(unique(Tcycles.Component));
    else
        S.NumComponentsInCycleTable = NaN;
    end

    % Edge-flow volumes.
    S.V_W = full(sum(Wwork(:)));
    S.V_C = full(sum(Cwork(:)));

    S.VolumeFromCycles = full(sum(lambda .* lengths));
    S.VolumeGap = S.VolumeFromCycles - S.V_C;
    S.VolumeAbsError = abs(S.VolumeGap);
    S.VolumeRelError = safeRelativeError(S.VolumeGap, S.V_C);

    % Completed circular throughput and mean cycle length.
    S.T_C = full(sum(lambda));

    if S.T_C > tol
        S.Lbar_C = S.V_C / S.T_C;
    else
        S.Lbar_C = NaN;
    end

    % Directional throughput from the original network imbalance.
    b = full(sum(Wwork, 2) - sum(Wwork, 1).');
    b(abs(b) <= tol) = 0;

    S.ImbalanceVector = b;
    S.ImbalanceResidual = full(abs(sum(b)));

    S.D_net = full(sum(max(b, 0)));
    S.D_net_abs_form = full(0.5 * sum(abs(b)));
    S.D_net_gap = S.D_net - S.D_net_abs_form;
    S.D_net_abs_error = abs(S.D_net_gap);

    % Circularity measures.
    if S.V_W > tol
        S.Gamma_vol = S.V_C / S.V_W;
    else
        S.Gamma_vol = NaN;
    end

    if S.T_C > tol && S.D_net > tol
        S.Gamma_thr = S.T_C / (S.T_C + S.D_net);
        S.R_thr = S.T_C / S.D_net;
    elseif S.T_C > tol && S.D_net <= tol
        S.Gamma_thr = 1;
        S.R_thr = Inf;
    elseif S.T_C <= tol && S.D_net > tol
        S.Gamma_thr = 0;
        S.R_thr = 0;
    else
        S.Gamma_thr = NaN;
        S.R_thr = NaN;
    end

    % Cycle-length summaries.
    if S.NumCycles > 0
        S.MinLength = min(lengths);
        S.MaxLength = max(lengths);
        S.UnweightedMeanLength = mean(lengths);
    else
        S.MinLength = NaN;
        S.MaxLength = NaN;
        S.UnweightedMeanLength = NaN;
    end

    S.ThroughputWeightedMeanLength = S.Lbar_C;

    % Cycle concentration.
    if S.T_C > tol
        p = lambda / S.T_C;
        S.HHI = sum(p.^2);
        S.N_eff = 1 / S.HHI;
    else
        p = zeros(size(lambda));
        S.HHI = NaN;
        S.N_eff = NaN;
    end

    % Top-k coverage tables.
    S.TopKByThroughput = makeTopKCoverage(lambda, lengths, S.T_C, S.V_C, ...
        opts.TopK, 'throughput');

    S.TopKByVolume = makeTopKCoverage(lambda, lengths, S.T_C, S.V_C, ...
        opts.TopK, 'volume');

    % Diagnostics.
    diagnostics = struct();

    diagnostics.BalanceResidualC = full(max(abs(sum(Cwork, 2) - sum(Cwork, 1).')));
    diagnostics.BalanceResidualW = full(max(abs(sum(Wwork, 2) - sum(Wwork, 1).')));
    diagnostics.ImbalanceSumResidual = S.ImbalanceResidual;

    diagnostics.VolumeIdentityChecked = opts.CheckVolumeIdentity;
    diagnostics.VolumeIdentityAbsError = S.VolumeAbsError;
    diagnostics.VolumeIdentityRelError = S.VolumeRelError;
    diagnostics.VolumeIdentityPass = true;

    if opts.CheckVolumeIdentity
        scale = max(1, S.V_C);
        diagnostics.VolumeIdentityPass = S.VolumeAbsError <= 1e-8 * scale;

        if ~diagnostics.VolumeIdentityPass && opts.ErrorOnCheckFailure
            error(['Cycle volume identity failed. abs(sum(lambda .* Length) - sum(C(:))) = %.3g. ' ...
                   'This may indicate incomplete enumeration or an estimated/partial cycle table.'], ...
                   S.VolumeAbsError);
        end
    end

    diagnostics.ReconstructionChecked = opts.CheckReconstruction;
    diagnostics.ReconstructionError = NaN;
    diagnostics.ReconstructionRelError = NaN;
    diagnostics.ReconstructionPass = true;

    if opts.CheckReconstruction
        Chat = reconstructFromCycleTable(n, Tcycles, opts);
        diagnostics.ReconstructionError = norm(Chat - Cwork, 'fro');
        diagnostics.ReconstructionRelError = diagnostics.ReconstructionError / max(1, norm(Cwork, 'fro'));
        diagnostics.ReconstructionPass = diagnostics.ReconstructionError <= 1e-8 * max(1, norm(Cwork, 'fro'));

        if ~diagnostics.ReconstructionPass && opts.ErrorOnCheckFailure
            error(['Cycle reconstruction failed. Frobenius error = %.3g. ' ...
                   'This may indicate incomplete enumeration, a partial Monte Carlo cycle table, ' ...
                   'or an inconsistent C/Tcycles pair.'], ...
                   diagnostics.ReconstructionError);
        end
    end

    diagnostics.DirectionalThroughputIdentityError = S.D_net_abs_error;
    diagnostics.DirectionalThroughputIdentityPass = S.D_net_abs_error <= 1e-8 * max(1, S.D_net_abs_form);

    if ~diagnostics.DirectionalThroughputIdentityPass && opts.ErrorOnCheckFailure
        error(['Directional throughput identity failed. abs(D_net - 0.5*sum(abs(b))) = %.3g.'], ...
               S.D_net_abs_error);
    end
end


function coverageTable = makeTopKCoverage(lambda, lengths, T_C, V_C, topK, rankingMode)
% makeTopKCoverage
% Top-k throughput and volume coverage under a chosen ranking.

    numCycles = numel(lambda);

    if numCycles == 0
        coverageTable = table();
        return;
    end

    volumeWeights = lambda(:) .* lengths(:);

    switch lower(rankingMode)
        case 'throughput'
            [~, order] = sort(lambda, 'descend');
        case 'volume'
            [~, order] = sort(volumeWeights, 'descend');
        otherwise
            error('Unknown rankingMode: %s', rankingMode);
    end

    topK = unique(topK(:));
    topK = topK(topK >= 1);
    topK = min(topK, numCycles);
    topK = unique(topK);

    K = zeros(numel(topK), 1);
    ThroughputCoverage = NaN(numel(topK), 1);
    VolumeCoverage = NaN(numel(topK), 1);

    cumLambda = cumsum(lambda(order));
    cumVolume = cumsum(volumeWeights(order));

    for r = 1:numel(topK)
        k = topK(r);
        K(r) = k;

        if T_C > 0
            ThroughputCoverage(r) = cumLambda(k) / T_C;
        end

        if V_C > 0
            VolumeCoverage(r) = cumVolume(k) / V_C;
        end
    end

    coverageTable = table(K, ThroughputCoverage, VolumeCoverage);
end


function Chat = reconstructFromCycleTable(n, Tcycles, opts)
% reconstructFromCycleTable
% Reconstruct circular component from cycle table.

    Chat = zeros(n);

    if isempty(Tcycles)
        return;
    end

    cycles = Tcycles.(opts.CycleVariable);
    lambda = full(Tcycles.(opts.WeightVariable));

    for r = 1:height(Tcycles)
        cyc = cycles{r};
        lam = lambda(r);

        if isempty(cyc)
            continue;
        end

        cyc = cyc(:).';
        closed = [cyc cyc(1)];

        for k = 1:numel(cyc)
            i = closed(k);
            j = closed(k+1);

            if i < 1 || i > n || j < 1 || j > n
                error('Cycle table contains node index outside the dimensions of W and C.');
            end

            Chat(i, j) = Chat(i, j) + lam;
        end
    end
end


function r = safeRelativeError(gap, scaleQuantity)
% safeRelativeError
% Relative error with robust denominator.

    r = abs(gap) / max(1, abs(scaleQuantity));
end


function validateSquareFinite(X, name)
% validateSquareFinite
% Basic square finite matrix validation.

    if ~ismatrix(X) || size(X, 1) ~= size(X, 2)
        error('%s must be a square matrix.', name);
    end

    if any(~isfinite(X(:)))
        error('%s contains non-finite entries.', name);
    end
end


function validateCycleTable(Tcycles, opts)
% validateCycleTable
% Validate required cycle-table variables.

    if ~istable(Tcycles)
        error('Tcycles must be a MATLAB table.');
    end

    required = {opts.CycleVariable, opts.LengthVariable, opts.WeightVariable};

    for k = 1:numel(required)
        if ~ismember(required{k}, Tcycles.Properties.VariableNames)
            error('Tcycles is missing required variable: %s', required{k});
        end
    end

    cycles = Tcycles.(opts.CycleVariable);

    if ~iscell(cycles)
        error('Tcycles.%s must be a cell array of node-index vectors.', opts.CycleVariable);
    end
end


function opts = parseCycleSummaryOptions(varargin)
% parseCycleSummaryOptions
% Parse name-value options.

    opts = struct();
    opts.Tol = 1e-12;
    opts.CheckVolumeIdentity = true;
    opts.CheckReconstruction = true;
    opts.ErrorOnCheckFailure = false;
    opts.TopK = [1 5 10 25 50 100];
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
            case 'tol'
                opts.Tol = value;

            case 'checkvolumeidentity'
                opts.CheckVolumeIdentity = value;

            case 'checkreconstruction'
                opts.CheckReconstruction = value;

            case 'erroroncheckfailure'
                opts.ErrorOnCheckFailure = value;

            case 'topk'
                opts.TopK = value;

            case 'weightvariable'
                opts.WeightVariable = value;

            case 'lengthvariable'
                opts.LengthVariable = value;

            case 'cyclevariable'
                opts.CycleVariable = value;

            otherwise
                error('Unknown option: %s', name);
        end
    end
end