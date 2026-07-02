function signed = cdfdSignedCycleSummaries(Tcycles, Sigma, varargin)
% cdfdSignedCycleSummaries
% Signed feedback summaries from probabilistic cycle coefficients.
%
% PURPOSE
%   Given a cycle table with probabilistic cycle coefficients lambda_gamma
%   and a sign matrix Sigma, compute the signed feedback decomposition:
%
%       T_C = T^+ + T^-
%       V_C = V^+ + V^-
%
%   where a cycle is positive if the product of signs around the cycle is
%   +1, and negative if the product is -1.
%
% INPUT
%   Tcycles   table with at least columns:
%               Cycle, Length, lambda
%
%             Usually returned by cdfdExactAllCycles(C).
%
%   Sigma     n x n sign matrix. Supported cycle edges must have entries
%             +1 or -1. Entries outside the cycle support may be zero.
%
% OPTIONS
%   'Tol'              default 1e-12
%   'TargetEdges'      k x 2 matrix of directed edges [i j] for contextual
%                      edge-role summaries. Default [].
%   'WeightVariable'   default 'lambda'
%   'LengthVariable'   default 'Length'
%   'CycleVariable'    default 'Cycle'
%
% OUTPUT
%   signed    struct with fields:
%
%             Cycle-level outputs:
%               CycleTable          original cycle table with added columns:
%                                     Polarity
%                                     IsPositive
%                                     IsNegative
%                                     ThroughputContribution
%                                     VolumeContribution
%                                     SignedThroughputContribution
%                                     SignedVolumeContribution
%
%             Aggregate throughput:
%               Tplus               sum lambda over positive cycles
%               Tminus              sum lambda over negative cycles
%               Ttotal              Tplus + Tminus
%               ShareTplus          Tplus / Ttotal
%               ShareTminus         Tminus / Ttotal
%
%             Aggregate edge-volume:
%               Vplus               sum lambda*Length over positive cycles
%               Vminus              sum lambda*Length over negative cycles
%               Vtotal              Vplus + Vminus
%               ShareVplus          Vplus / Vtotal
%               ShareVminus         Vminus / Vtotal
%
%             Mean signed cycle lengths:
%               LbarPlus            Vplus / Tplus
%               LbarMinus           Vminus / Tminus
%               LbarTotal           Vtotal / Ttotal
%
%             Counts:
%               NumCycles
%               NumPositiveCycles
%               NumNegativeCycles
%
%             Edge-role outputs:
%               Cplus               matrix with c_e^+ values
%               Cminus              matrix with c_e^- values
%               Ctotal              Cplus + Cminus
%               EdgeRoleTable       table with one row per edge appearing in
%                                    at least one listed cycle
%
%             Target-edge outputs:
%               TargetEdges         requested target edges
%               TargetEdgeTable     contextual role table for target edges
%
%             Convenience fields, when exactly one target edge is supplied:
%               TargetEdge
%               TargetEdgeCplus
%               TargetEdgeCminus
%               TargetEdgeCtotal
%
% NOTES
%   - Self-loops are allowed and are treated as length-one cycles.
%   - A self-loop at node i has polarity Sigma(i,i).
%   - The function assumes Tcycles is a complete or intended cycle table.
%     If Tcycles is a partial Monte Carlo table, the summaries are summaries
%     of the listed/estimated cycles only.

    opts = parseSignedOptions(varargin{:});
    tol = opts.Tol;

    validateCycleTable(Tcycles, opts);
    validateSigma(Sigma, tol);

    Sigma = full(double(Sigma));
    Sigma(abs(Sigma) <= tol) = 0;
    Sigma(Sigma > 0) = 1;
    Sigma(Sigma < 0) = -1;

    n = size(Sigma, 1);

    lambda = full(Tcycles.(opts.WeightVariable));
    lengths = full(Tcycles.(opts.LengthVariable));
    cycles = Tcycles.(opts.CycleVariable);

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

    numCycles = height(Tcycles);

    polarity = zeros(numCycles, 1);
    volumeContribution = lambda .* lengths;

    Cplus = zeros(n);
    Cminus = zeros(n);

    for r = 1:numCycles
        cyc = normaliseCycle(cycles{r});

        if isempty(cyc)
            error('Cycle table contains an empty cycle.');
        end

        if any(cyc < 1) || any(cyc > n) || any(cyc ~= round(cyc))
            error('Cycle table contains node index outside the dimensions of Sigma.');
        end

        if numel(unique(cyc)) ~= numel(cyc)
            error('Cycle table contains a non-simple cycle or repeated node.');
        end

        if lengths(r) ~= numel(cyc)
            error('Cycle length mismatch in row %d.', r);
        end

        closed = [cyc cyc(1)];

        sig = 1;

        for k = 1:numel(cyc)
            i = closed(k);
            j = closed(k+1);

            edgeSign = Sigma(i, j);

            if edgeSign == 0
                error('Sigma(%d,%d) is zero, but edge %d -> %d appears in cycle row %d.', ...
                    i, j, i, j, r);
            end

            sig = sig * edgeSign;
        end

        polarity(r) = sig;

        for k = 1:numel(cyc)
            i = closed(k);
            j = closed(k+1);

            if sig == 1
                Cplus(i, j) = Cplus(i, j) + lambda(r);
            elseif sig == -1
                Cminus(i, j) = Cminus(i, j) + lambda(r);
            else
                error('Unexpected cycle polarity in row %d.', r);
            end
        end
    end

    isPositive = polarity == 1;
    isNegative = polarity == -1;

    Tsigned = Tcycles;
    Tsigned.Polarity = polarity;
    Tsigned.IsPositive = isPositive;
    Tsigned.IsNegative = isNegative;
    Tsigned.ThroughputContribution = lambda;
    Tsigned.VolumeContribution = volumeContribution;
    Tsigned.SignedThroughputContribution = polarity .* lambda;
    Tsigned.SignedVolumeContribution = polarity .* volumeContribution;

    signed = struct();

    signed.CycleTable = Tsigned;

    signed.Tplus = full(sum(lambda(isPositive)));
    signed.Tminus = full(sum(lambda(isNegative)));
    signed.Ttotal = signed.Tplus + signed.Tminus;

    signed.Vplus = full(sum(volumeContribution(isPositive)));
    signed.Vminus = full(sum(volumeContribution(isNegative)));
    signed.Vtotal = signed.Vplus + signed.Vminus;

    signed.ShareTplus = safeShare(signed.Tplus, signed.Ttotal, tol);
    signed.ShareTminus = safeShare(signed.Tminus, signed.Ttotal, tol);

    signed.ShareVplus = safeShare(signed.Vplus, signed.Vtotal, tol);
    signed.ShareVminus = safeShare(signed.Vminus, signed.Vtotal, tol);

    signed.LbarPlus = safeRatio(signed.Vplus, signed.Tplus, tol);
    signed.LbarMinus = safeRatio(signed.Vminus, signed.Tminus, tol);
    signed.LbarTotal = safeRatio(signed.Vtotal, signed.Ttotal, tol);

    signed.NumCycles = numCycles;
    signed.NumPositiveCycles = full(sum(isPositive));
    signed.NumNegativeCycles = full(sum(isNegative));

    Cplus(abs(Cplus) <= tol) = 0;
    Cminus(abs(Cminus) <= tol) = 0;

    signed.Cplus = Cplus;
    signed.Cminus = Cminus;
    signed.Ctotal = Cplus + Cminus;

    signed.EdgeRoleTable = makeEdgeRoleTable(Cplus, Cminus, tol);

    targetEdges = opts.TargetEdges;

    if isempty(targetEdges)
        signed.TargetEdges = zeros(0, 2);
        signed.TargetEdgeTable = table();
    else
        validateTargetEdges(targetEdges, n);
        signed.TargetEdges = targetEdges;
        signed.TargetEdgeTable = makeTargetEdgeTable(Cplus, Cminus, targetEdges, tol);

        if size(targetEdges, 1) == 1
            i = targetEdges(1, 1);
            j = targetEdges(1, 2);

            signed.TargetEdge = targetEdges;
            signed.TargetEdgeCplus = Cplus(i, j);
            signed.TargetEdgeCminus = Cminus(i, j);
            signed.TargetEdgeCtotal = Cplus(i, j) + Cminus(i, j);
        end
    end

    signed.Options = opts;
end


function edgeTable = makeEdgeRoleTable(Cplus, Cminus, tol)
% makeEdgeRoleTable
% Build table of contextual positive/negative edge roles.

    Ctotal = Cplus + Cminus;

    [from, to, ctotal] = find(Ctotal);

    if isempty(from)
        edgeTable = table();
        return;
    end

    numEdges = numel(from);

    cplus = zeros(numEdges, 1);
    cminus = zeros(numEdges, 1);

    for r = 1:numEdges
        cplus(r) = Cplus(from(r), to(r));
        cminus(r) = Cminus(from(r), to(r));
    end

    positiveShare = NaN(numEdges, 1);
    negativeShare = NaN(numEdges, 1);

    positiveShare(ctotal > tol) = cplus(ctotal > tol) ./ ctotal(ctotal > tol);
    negativeShare(ctotal > tol) = cminus(ctotal > tol) ./ ctotal(ctotal > tol);

    edgeTable = table( ...
        from, ...
        to, ...
        cplus, ...
        cminus, ...
        ctotal, ...
        positiveShare, ...
        negativeShare, ...
        'VariableNames', { ...
            'From', ...
            'To', ...
            'Cplus', ...
            'Cminus', ...
            'Ctotal', ...
            'PositiveShare', ...
            'NegativeShare' ...
        } ...
    );

    edgeTable = sortrows(edgeTable, {'From', 'To'});
end


function targetTable = makeTargetEdgeTable(Cplus, Cminus, targetEdges, tol)
% makeTargetEdgeTable
% Build contextual role table for requested target edges.

    kEdges = size(targetEdges, 1);

    from = targetEdges(:, 1);
    to = targetEdges(:, 2);

    cplus = zeros(kEdges, 1);
    cminus = zeros(kEdges, 1);
    ctotal = zeros(kEdges, 1);

    for r = 1:kEdges
        i = from(r);
        j = to(r);

        cplus(r) = Cplus(i, j);
        cminus(r) = Cminus(i, j);
        ctotal(r) = cplus(r) + cminus(r);
    end

    positiveShare = NaN(kEdges, 1);
    negativeShare = NaN(kEdges, 1);

    positiveShare(ctotal > tol) = cplus(ctotal > tol) ./ ctotal(ctotal > tol);
    negativeShare(ctotal > tol) = cminus(ctotal > tol) ./ ctotal(ctotal > tol);

    targetTable = table( ...
        from, ...
        to, ...
        cplus, ...
        cminus, ...
        ctotal, ...
        positiveShare, ...
        negativeShare, ...
        'VariableNames', { ...
            'From', ...
            'To', ...
            'Cplus', ...
            'Cminus', ...
            'Ctotal', ...
            'PositiveShare', ...
            'NegativeShare' ...
        } ...
    );
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


function y = safeShare(x, total, tol)

    if total > tol
        y = x / total;
    else
        y = NaN;
    end
end


function y = safeRatio(num, den, tol)

    if den > tol
        y = num / den;
    else
        y = NaN;
    end
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


function validateSigma(Sigma, tol)

    if ~ismatrix(Sigma) || size(Sigma, 1) ~= size(Sigma, 2)
        error('Sigma must be a square sign matrix.');
    end

    if any(~isfinite(Sigma(:)))
        error('Sigma contains non-finite entries.');
    end

    Sigma = full(double(Sigma));
    Sigma(abs(Sigma) <= tol) = 0;

    nonzero = Sigma ~= 0;
    invalid = nonzero & abs(abs(Sigma) - 1) > tol;

    if any(invalid(:))
        error('Non-zero entries of Sigma must be +1 or -1.');
    end
end


function validateTargetEdges(targetEdges, n)

    if ~isnumeric(targetEdges) || size(targetEdges, 2) ~= 2
        error('TargetEdges must be a k x 2 numeric matrix of directed edges [i j].');
    end

    if any(targetEdges(:) < 1) || any(targetEdges(:) > n) || any(targetEdges(:) ~= round(targetEdges(:)))
        error('TargetEdges contains invalid node indices.');
    end
end


function opts = parseSignedOptions(varargin)

    opts = struct();
    opts.Tol = 1e-12;
    opts.TargetEdges = zeros(0, 2);
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

            case 'targetedges'
                opts.TargetEdges = value;

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