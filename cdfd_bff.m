function [C, D, info] = cdfd_bff(W, varargin)
%CDFD_BFF  Balanced Flow Forwarding (BFF) CDFD decomposition.
%
%   [C, D, info] = CDFD_BFF(W) decomposes a weighted directed adjacency
%   matrix W into:
%       C  - circular component
%       D  - directional (acyclic) component
%
%   This is a MATLAB-first implementation of the recursive BFF idea:
%     1. remove self-loops (which belong in C),
%     2. remove SCC-singletons from the active subproblem,
%     3. on each remaining nontrivial SCC, compute one BFF circulation step,
%     4. subtract and repeat until the residual is acyclic.
%
%   Name-value options:
%       'ToleranceZero'  : absolute zero tolerance (default 1e-12)
%       'MaxIterations'  : max recursive iterations (default 2*nnz(W)+1)
%       'Validate'       : true/false, run decomposition checks (default true)
%
%   Notes
%   -----
%   - W must be square, finite, and nonnegative.
%   - Rows are sources, columns are destinations.
%   - Self-loops are assigned to C automatically.
%
%   Output info fields:
%       .iterations
%       .maxIterations
%       .balanceError
%       .sumError
%       .isAcyclic
%       .circularity
%       .directionality
%
%   See also: DIGRAPH, CONNCOMP, EIGS

    p = inputParser;
    addParameter(p, 'ToleranceZero', 1e-12, @(x) isnumeric(x) && isscalar(x) && x > 0);
    addParameter(p, 'MaxIterations', [], @(x) isempty(x) || (isnumeric(x) && isscalar(x) && x >= 0));
    addParameter(p, 'Validate', true, @(x) islogical(x) && isscalar(x));
    parse(p, varargin{:});

    tol = p.Results.ToleranceZero;
    validateOut = p.Results.Validate;

    % ---- checks ----
    if ~ismatrix(W) || size(W,1) ~= size(W,2)
        error('cdfd_bff:WNotSquare', 'W must be a square adjacency matrix.');
    end
    if any(~isfinite(W(:)))
        error('cdfd_bff:NonFinite', 'W must contain only finite values.');
    end
    if any(W(:) < -tol)
        error('cdfd_bff:NegativeWeights', 'W must be nonnegative.');
    end

    W = sparse(double(W));
    W = pruneSmall_(W, tol);

    n = size(W,1);

    % ---- self-loops belong to C ----
    loops = spdiags(diag(W), 0, n, n);
    Wnol = W - loops;
    Wnol = pruneSmall_(Wnol, tol);

    % ---- separate SCC-singletons from active problem ----
    [Wcore, Wisolated, idxCore] = separateIsolatedSCCSingletons_(Wnol, tol);

    Wtemp = Wcore;
    nCore = size(Wtemp,1);
    nEdges = nnz(Wtemp);

    if isempty(p.Results.MaxIterations)
        maxIt = 2 * nEdges + 1;
    else
        maxIt = p.Results.MaxIterations;
    end

    iterations = 0;
    nComponents = 0;

    while nComponents < nCore && iterations < maxIt
        [cStep, nComponents] = bffStep_(Wtemp, tol);
        Wtemp = Wtemp - cStep;
        Wtemp = pruneSmall_(Wtemp, tol);
        iterations = iterations + 1;
    end

    if nComponents < nCore
        error('cdfd_bff:MaxIterationsReached', ...
            'Maximum iterations reached before obtaining an acyclic residual.');
    end

    % ---- rebuild full residual D ----
    D = sparse(n, n);
    if ~isempty(idxCore)
        D(idxCore, idxCore) = Wtemp;
    end
    D = D + Wisolated;
    D = pruneSmall_(D, tol);

    % ---- recover C on full matrix (including loops) ----
    C = W - D;
    C = pruneSmall_(C, tol);

    % ---- diagnostics ----
    info = struct();
    info.iterations    = iterations;
    info.maxIterations = maxIt;

    checks = decompositionChecks_(C, D, W, tol);
    info.balanceError  = checks.balanceError;
    info.sumError      = checks.sumError;
    info.isAcyclic     = checks.isAcyclic;

    totalW = full(sum(W(:)));
    totalC = full(sum(C(:)));
    if totalW > 0
        info.circularity   = totalC / totalW;
        info.directionality = 1 - info.circularity;
    else
        info.circularity   = 0;
        info.directionality = 0;
    end

    if validateOut
        if ~checks.isAcyclic
            warning('cdfd_bff:NonAcyclicResidual', ...
                'Directional part is not acyclic.');
        end
        if checks.balanceError > 20 * nnz(W) * tol
            warning('cdfd_bff:ImbalancedCircularPart', ...
                'Circular part balance error is %g.', checks.balanceError);
        end
        if checks.sumError > 20 * nnz(W) * tol
            warning('cdfd_bff:DoesNotSumToW', ...
                'C + D differs from W by %g.', checks.sumError);
        end
    end
end

% =====================================================================
function [Cstep, nComponents] = bffStep_(W, tol)
%BFFSTEP_ One BFF step on current residual W.

    n = size(W,1);
    if n == 0
        Cstep = sparse(0,0);
        nComponents = 0;
        return;
    end

    G = digraph(spones(W));
    bins = conncomp(G, 'Type', 'strong');
    nComponents = max(bins);

    counts = accumarray(bins(:), 1, [nComponents 1]);
    nontrivial = find(counts > 1);   % loops were removed already

    Cstep = sparse(n, n);

    for k = 1:numel(nontrivial)
        compID = nontrivial(k);
        idx = find(bins == compID);
        Wsub = W(idx, idx);

        Csub = bffStronglyConnected_(Wsub, tol);
        Cstep(idx, idx) = Csub;
    end

    Cstep = pruneSmall_(Cstep, tol);
end

% =====================================================================
function Csub = bffStronglyConnected_(W, tol)
%BFFSTRONGLYCONNECTED_ One BFF circulation on a strongly connected block.

    n = size(W,1);
    if n == 0
        Csub = sparse(0,0);
        return;
    elseif n == 1
        % Should not happen here because loops removed and SCC-singletons excluded,
        % but keep harmless fallback.
        Csub = sparse(1,1);
        return;
    end

    out = full(sum(W, 2));
    if any(out <= tol)
        error('cdfd_bff:ZeroOutflowInSCC', ...
            'Encountered a strongly connected block with zero outflow node.');
    end

    P = spdiags(1 ./ out, 0, n, n) * W;   % row-stochastic
    pi = stationaryDistribution_(P, tol); % column vector, sums to 1

    Cstationary = spdiags(pi, 0, n, n) * P;

    mask = Cstationary > tol;
    if ~any(mask(:))
        Csub = sparse(n,n);
        return;
    end

    ratios = full(W(mask) ./ Cstationary(mask));
    alpha = min(ratios);

    Csub = alpha * Cstationary;
    Csub = pruneSmall_(Csub, tol);
end

% =====================================================================
function pi = stationaryDistribution_(P, tol)
%STATIONARYDISTRIBUTION_ Left stationary distribution of row-stochastic P.

    n = size(P,1);

    if n == 1
        pi = 1;
        return;
    end

    % Solve (P' - I) pi = 0 with sum(pi)=1 by replacing one row.
    A = P.' - speye(n);
    b = zeros(n,1);
    A(end,:) = 1;
    b(end) = 1;

    pi = A \ b;
    pi = real(pi);
    pi(pi < 0) = 0;

    s = sum(pi);
    if s <= tol
        % Fallback: simple power iteration on row vector
        x = ones(1,n) / n;
        for t = 1:5000
            xNew = x * P;
            if norm(xNew - x, 1) < 1e-13
                break;
            end
            x = xNew;
        end
        pi = x(:);
        pi(pi < 0) = 0;
        s = sum(pi);
        if s <= tol
            error('cdfd_bff:StationaryDistributionFailed', ...
                'Failed to compute a stationary distribution.');
        end
    end

    pi = pi / s;
end

% =====================================================================
function [Wcore, Wisolated, idxCore] = separateIsolatedSCCSingletons_(W, tol)
%SEPARATEISOLATEDSCCSINGLETONS_ Remove SCC-singletons from active recursion.

    n = size(W,1);
    if n == 0
        Wcore = sparse(0,0);
        Wisolated = sparse(0,0);
        idxCore = [];
        return;
    end

    G = digraph(spones(W));
    bins = conncomp(G, 'Type', 'strong');
    nComponents = max(bins);
    counts = accumarray(bins(:), 1, [nComponents 1]);

    keepComp = find(counts > 1);
    idxCore = find(ismember(bins, keepComp));

    Wcore = W(idxCore, idxCore);

    Wisolated = sparse(n, n);
    if ~isempty(idxCore)
        Wembedded = sparse(n, n);
        Wembedded(idxCore, idxCore) = Wcore;
        Wisolated = W - Wembedded;
    else
        Wisolated = W;
    end

    Wcore = pruneSmall_(Wcore, tol);
    Wisolated = pruneSmall_(Wisolated, tol);
end

% =====================================================================
function checks = decompositionChecks_(C, D, W, tol)
%DECOMPOSITIONCHECKS_ Basic validity diagnostics.

    checks = struct();

    checks.balanceError = norm(full(sum(C,2) - sum(C,1).'), 2);
    checks.sumError = norm(W - (C + D), 'fro');

    if any(diag(D) > tol)
        checks.isAcyclic = false;
        return;
    end

    if size(D,1) == 0
        checks.isAcyclic = true;
        return;
    end

    G = digraph(spones(D));
    bins = conncomp(G, 'Type', 'strong');
    counts = accumarray(bins(:), 1);
    checks.isAcyclic = all(counts == 1);
end

% =====================================================================
function X = pruneSmall_(X, tol)
%PRUNESMALL_ Zero out tiny entries in sparse matrix.
    X = sparse(X);
    if nnz(X) == 0
        return;
    end
    mask = abs(X) < tol;
    if any(mask)
        X(mask) = 0;
        X = sparse(X);
    end
end