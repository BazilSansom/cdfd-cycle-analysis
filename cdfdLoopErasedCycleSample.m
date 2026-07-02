function sample = cdfdLoopErasedCycleSample(A, varargin)
% cdfdLoopErasedCycleSample
% Generate one chronological loop-erased cycle sample from a Markov chain.
%
% PURPOSE
%   Simulates a Markov chain until the current loop-erased path closes a
%   simple directed cycle. This is the elementary sampling step behind the
%   Qian--Qian / Kalpazidou probabilistic cycle representation.
%
%   The function maintains a self-avoiding path. At each Markov transition,
%   if the next state is new, it is appended to the path. If the next state
%   already appears in the path, the segment from that previous occurrence
%   to the current end of the path is returned as a simple cycle, and the
%   path is updated by erasing that cycle.
%
% INPUT
%   A   Either:
%         - a row-stochastic transition matrix P, or
%         - a non-negative flow/circulation matrix C, from which
%           P_ij = C_ij / sum_j C_ij is constructed.
%
% OPTIONS
%   'InputType'          default 'auto'. One of:
%                          'auto'        infer transition vs flow
%                          'transition'  treat A as row-stochastic P
%                          'flow'        row-normalise A into P
%   'Tol'                default 1e-12
%   'StartNode'          default []. Starting node if InitialPath is empty.
%   'InitialPath'        default []. Self-avoiding current path. If supplied,
%                        the last node is treated as the current Markov state.
%   'MaxSteps'           default 1e6. Maximum Markov transitions before error.
%   'ReturnClosedCycle'  default true. Include closed representation
%                        [cycle cycle(1)] in sample.CycleClosed.
%   'RandomStream'       default []. Optional RandStream object.
%
% OUTPUT
%   sample   struct with fields:
%              Cycle              open cycle node list, e.g. [1 2 3]
%              CycleClosed        closed cycle, e.g. [1 2 3 1], or []
%              CanonicalCycle     cycle up to cyclic rotation
%              Key                string key for CanonicalCycle
%              Length             cycle length
%              IsSelfLoop         true for length-one self-loop cycles
%              InitialPath        input starting path
%              UpdatedPath        path after erasing the sampled cycle
%              StartNode          first node of InitialPath
%              EndNode            last node of UpdatedPath
%              TransitionCount    Markov transitions used to form the cycle
%              InputTypeUsed      'transition' or 'flow'
%              Completed          true if a cycle was sampled
%              Options            parsed options
%
% NOTES
%   - Self-loops are allowed. A transition i -> i when the current path
%     ends at i returns the length-one cycle [i].
%   - Cycles are returned in open representation. The canonical cycle is
%     unique up to cyclic rotation but does not reverse orientation.
%   - For repeated streaming use, pass sample.UpdatedPath as the next
%     InitialPath. This preserves the chronological loop-erasure state.

    opts = parseLoopErasedOptions(varargin{:});
    tol = opts.Tol;

    validateSquareNonnegative(A, tol);

    [P, inputTypeUsed, startProb] = makeTransitionMatrix(A, opts);

    n = size(P, 1);

    if isempty(opts.InitialPath)
        if isempty(opts.StartNode)
            startNode = sampleDiscrete(startProb, opts.RandomStream);
        else
            startNode = opts.StartNode;
            validateNode(startNode, n, 'StartNode');
        end

        path = startNode;
    else
        path = opts.InitialPath(:).';
        validatePath(path, n);

        if ~isempty(opts.StartNode) && opts.StartNode ~= path(1)
            error('If both StartNode and InitialPath are supplied, StartNode must equal InitialPath(1).');
        end

        startNode = path(1);
    end

    if numel(unique(path)) ~= numel(path)
        error('InitialPath must be self-avoiding, with no repeated nodes.');
    end

    initialPath = path;

    completed = false;
    cycle = [];
    updatedPath = path;
    transitionCount = 0;

    for step = 1:opts.MaxSteps
        current = path(end);

        next = sampleNextState(P, current, tol, opts.RandomStream);
        transitionCount = transitionCount + 1;

        previousPosition = find(path == next, 1, 'first');

        if isempty(previousPosition)
            path(end+1) = next; %#ok<AGROW>
        else
            cycle = path(previousPosition:end);
            updatedPath = path(1:previousPosition);
            completed = true;
            break;
        end
    end

    if ~completed
        error('No loop-erased cycle was formed within MaxSteps=%d.', opts.MaxSteps);
    end

    canonicalCycle = canonicalDirectedCycle(cycle);
    key = cycleKey(canonicalCycle);

    if opts.ReturnClosedCycle
        cycleClosed = [cycle cycle(1)];
    else
        cycleClosed = [];
    end

    sample = struct();
    sample.Cycle = cycle;
    sample.CycleClosed = cycleClosed;
    sample.CanonicalCycle = canonicalCycle;
    sample.Key = key;
    sample.Length = numel(cycle);
    sample.IsSelfLoop = numel(cycle) == 1;
    sample.InitialPath = initialPath;
    sample.UpdatedPath = updatedPath;
    sample.StartNode = startNode;
    sample.EndNode = updatedPath(end);
    sample.TransitionCount = transitionCount;
    sample.InputTypeUsed = inputTypeUsed;
    sample.Completed = completed;
    sample.Options = opts;
end


function [P, inputTypeUsed, startProb] = makeTransitionMatrix(A, opts)

    tol = opts.Tol;
    A = sparse(double(A));
    A(abs(A) <= tol) = 0;
    A = sparse(A);

    n = size(A, 1);
    rowSum = full(sum(A, 2));

    inputType = lower(string(opts.InputType));

    if inputType == "auto"
        positiveRows = rowSum > tol;

        if any(positiveRows) && max(abs(rowSum(positiveRows) - 1)) <= 1e-8
            inputType = "transition";
        else
            inputType = "flow";
        end
    end

    switch inputType
        case "transition"
            positiveRows = rowSum > tol;

            if ~all(positiveRows)
                error('Transition matrix has one or more zero-outflow rows.');
            end

            if max(abs(rowSum - 1)) > 1e-8
                error('InputType is transition, but row sums are not one within tolerance.');
            end

            P = spdiags(1 ./ rowSum, 0, n, n) * A;
            inputTypeUsed = 'transition';

            % Without stationary information, use a uniform start over all nodes.
            startProb = ones(n, 1) / n;

        case "flow"
            positiveRows = rowSum > tol;

            if ~any(positiveRows)
                error('Flow matrix has zero total outflow.');
            end

            invRow = zeros(n, 1);
            invRow(positiveRows) = 1 ./ rowSum(positiveRows);

            P = spdiags(invRow, 0, n, n) * A;
            inputTypeUsed = 'flow';

            % For a balanced circulation this is the stationary distribution.
            startProb = rowSum / sum(rowSum);

        otherwise
            error('InputType must be auto, transition, or flow.');
    end
end


function next = sampleNextState(P, current, tol, stream)

    probs = full(P(current, :));
    probs(abs(probs) <= tol) = 0;

    s = sum(probs);

    if s <= tol
        error('Current node %d has zero outgoing transition probability.', current);
    end

    probs = probs / s;

    u = randUniform(stream);

    cs = cumsum(probs);
    next = find(u <= cs, 1, 'first');

    if isempty(next)
        % Numerical fallback when u is extremely close to one.
        next = find(probs > 0, 1, 'last');
    end
end


function idx = sampleDiscrete(prob, stream)

    prob = full(prob(:));
    prob(prob < 0) = 0;

    s = sum(prob);

    if s <= 0
        error('Cannot sample from a zero probability vector.');
    end

    prob = prob / s;

    u = randUniform(stream);
    cs = cumsum(prob);

    idx = find(u <= cs, 1, 'first');

    if isempty(idx)
        idx = find(prob > 0, 1, 'last');
    end
end


function u = randUniform(stream)

    if isempty(stream)
        u = rand();
    else
        u = rand(stream);
    end
end


function canonical = canonicalDirectedCycle(cyc)
% canonicalDirectedCycle
% Canonicalise a directed cycle up to cyclic rotation.
%
% The orientation is preserved. For example:
%   [2 3 1] -> [1 2 3]
%   [3 2 1] -> [1 3 2]
%
% These are different directed cycles.

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


function validateSquareNonnegative(A, tol)

    if ~ismatrix(A) || size(A, 1) ~= size(A, 2)
        error('Input matrix must be square.');
    end

    if any(~isfinite(A(:)))
        error('Input matrix contains non-finite entries.');
    end

    if any(A(:) < -tol)
        error('Input matrix must be non-negative up to numerical tolerance.');
    end
end


function validateNode(node, n, name)

    if ~isnumeric(node) || ~isscalar(node) || ~isfinite(node) || ...
            node < 1 || node > n || node ~= round(node)
        error('%s must be a valid node index.', name);
    end
end


function validatePath(path, n)

    if isempty(path)
        error('InitialPath must be non-empty if supplied.');
    end

    if any(path < 1) || any(path > n) || any(path ~= round(path))
        error('InitialPath contains invalid node indices.');
    end
end


function opts = parseLoopErasedOptions(varargin)

    opts = struct();
    opts.InputType = 'auto';
    opts.Tol = 1e-12;
    opts.StartNode = [];
    opts.InitialPath = [];
    opts.MaxSteps = 1e6;
    opts.ReturnClosedCycle = true;
    opts.RandomStream = [];

    if mod(numel(varargin), 2) ~= 0
        error('Options must be name-value pairs.');
    end

    for k = 1:2:numel(varargin)
        name = varargin{k};
        value = varargin{k+1};

        switch lower(name)
            case 'inputtype'
                opts.InputType = value;

            case 'tol'
                opts.Tol = value;

            case 'startnode'
                opts.StartNode = value;

            case 'initialpath'
                opts.InitialPath = value;

            case 'maxsteps'
                opts.MaxSteps = value;

            case 'returnclosedcycle'
                opts.ReturnClosedCycle = value;

            case 'randomstream'
                opts.RandomStream = value;

            otherwise
                error('Unknown option: %s.', name);
        end
    end
end