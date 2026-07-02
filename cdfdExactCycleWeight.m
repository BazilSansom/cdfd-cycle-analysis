function out = cdfdExactCycleWeight(C, cycleNodes, varargin)
% cdfdExactCycleWeight
% Exact Qian--Qian / Kalpazidou probabilistic cycle coefficient for a
% nominated simple directed cycle.
%
% INPUTS
%   C           n x n non-negative circular component, typically C_BFF.
%   cycleNodes  vector of node indices, e.g. [1], [1 2], [1 2 3],
%               or optionally closed as [1 1], [1 2 3 1].
%
% NAME-VALUE OPTIONS
%   'Tol'           numerical tolerance, default 1e-12.
%   'CheckBalance'  whether to warn if C is not balanced, default true.
%
% OUTPUT
%   out is a struct with fields:
%     .cycle          open cycle node list, e.g. [1 2 3], or [1] for a self-loop
%     .componentNodes nodes in the strongly connected component
%     .cycleLocal     local node indices inside the component
%     .length         cycle length
%     .isSelfLoop     true for length-one self-loop cycles
%     .S              circular edge-flow volume of the component
%     .q              chronological loop-erased cycle frequency
%     .lambda         cycle coefficient lambda = S*q
%     .logq           log(q), useful for diagnostics
%
% NOTES
%   - Self-loops are allowed and are treated as length-one simple cycles.
%   - For a self-loop at node i, pass cycleNodes = [i] or [i i].
%   - C must be non-negative. For signed causal networks, pass the BFF
%     circular component of the non-negative magnitude network.

    opts = parseOptions(varargin{:});
    tol = opts.Tol;

    validateSquareNonnegative(C, tol);

    % Store sparse internally, but return ordinary full scalar summaries.
    C = sparse(double(C));
    C(abs(C) <= tol) = 0;
    C = sparse(C);

    n = size(C, 1);

    cycleNodes = cycleNodes(:).';

    % Allow optional closure, including [i i] for a self-loop.
    if numel(cycleNodes) >= 2 && cycleNodes(end) == cycleNodes(1)
        cycleNodes = cycleNodes(1:end-1);
    end

    if isempty(cycleNodes)
        error('cycleNodes must contain at least one node.');
    end

    if any(cycleNodes < 1) || any(cycleNodes > n) || any(cycleNodes ~= round(cycleNodes))
        error('cycleNodes must be valid integer node indices.');
    end

    if numel(unique(cycleNodes)) ~= numel(cycleNodes)
        error('cycleNodes must be a simple cycle with no repeated nodes, except optional closure.');
    end

    if opts.CheckBalance
        balanceResidual = full(max(abs(sum(C, 2) - sum(C, 1).')));
        scale = full(max(1, max(sum(C, 2))));

        if balanceResidual > 100 * tol * scale
            warning('C does not appear balanced. Max balance residual = %.3g.', balanceResidual);
        end
    end

    % Check nominated cycle edges exist.
    closedCycle = [cycleNodes cycleNodes(1)];

    for k = 1:numel(cycleNodes)
        i = closedCycle(k);
        j = closedCycle(k+1);

        if full(C(i, j)) <= tol
            error('Cycle edge %d -> %d has zero or negligible weight in C.', i, j);
        end
    end

    % Strong components of positive support. Self-loops are included.
    G = digraph(C > tol);
    comp = conncomp(G, 'Type', 'strong');

    c0 = comp(cycleNodes(1));

    if any(comp(cycleNodes) ~= c0)
        error('All cycle nodes must lie in the same strongly connected component.');
    end

    componentNodes = find(comp == c0);
    Csub = C(componentNodes, componentNodes);

    % Map global node labels to local component labels.
    localMap = zeros(1, n);
    localMap(componentNodes) = 1:numel(componentNodes);
    cycleLocal = localMap(cycleNodes);

    [q, lambda, logq, S] = qqWeightLocal(Csub, cycleLocal, tol);

    out = struct();
    out.cycle = cycleNodes;
    out.componentNodes = componentNodes;
    out.cycleLocal = cycleLocal;
    out.length = numel(cycleNodes);
    out.isSelfLoop = numel(cycleNodes) == 1;
    out.S = S;
    out.q = q;
    out.lambda = lambda;
    out.logq = logq;
end


function [q, lambda, logq, S] = qqWeightLocal(Csub, cycleLocal, tol)
% qqWeightLocal
% Compute q_gamma and lambda_gamma for a cycle inside one SCC.
%
% cycleLocal is an open row vector of local node indices.
% A self-loop is represented by a length-one vector, e.g. [i].

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

    logq = full(logProdP + logMinorNumerator - logDenominator);
    q = full(exp(logq));
    lambda = full(S * q);
end


function ld = logdetAbs(A)
% logdetAbs
% Log absolute determinant, with det(empty)=1.
% For the killed-chain minors used here, determinants should be positive.
% We use log(abs(det)) for numerical stability.

    if isempty(A)
        ld = 0;
        return;
    end

    A = full(A);

    [~, U, ~] = lu(A, 'vector');
    diagU = diag(U);

    if any(abs(diagU) == 0)
        ld = -inf;
    else
        ld = full(sum(log(abs(diagU))));
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

    y = full(y);
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


function opts = parseOptions(varargin)

    opts = struct();
    opts.Tol = 1e-12;
    opts.CheckBalance = true;

    if mod(numel(varargin), 2) ~= 0
        error('Options must be name-value pairs.');
    end

    for k = 1:2:numel(varargin)
        name = varargin{k};
        value = varargin{k+1};

        switch lower(name)
            case 'tol'
                opts.Tol = value;

            case 'checkbalance'
                opts.CheckBalance = value;

            otherwise
                error('Unknown option: %s', name);
        end
    end
end