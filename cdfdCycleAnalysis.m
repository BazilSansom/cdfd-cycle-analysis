function R = cdfdCycleAnalysis(W, varargin)
% cdfdCycleAnalysis
% Full exact CDFD cycle analysis pipeline.
%
% PURPOSE
%   Runs BFF CDFD, exact probabilistic all-cycle decomposition of the BFF
%   circular component, and network-level throughput/circularity summaries.
%
% WORKFLOW
%       W
%       -> cdfd_bff(W)
%       -> C, D
%       -> cdfdExactAllCycles(C)
%       -> Tcycles
%       -> cdfdCycleSummaries(W, C, Tcycles)
%       -> summary statistics
%
% INPUT
%   W   n x n non-negative weighted directed adjacency matrix.
%
% OPTIONS
%   'Tol'                  default 1e-12
%   'MaxCycles'            default 50000
%   'CheckBalance'         default true
%   'CheckReconstruction'  default true
%   'BFFValidate'          default true
%
% OUTPUT
%   R   struct with fields:
%         W              original network
%         C              BFF circular component
%         D              BFF directional component
%         bffInfo        diagnostics from cdfd_bff
%         Tcycles        exact probabilistic cycle table
%         cycleStatus    diagnostics from cdfdExactAllCycles
%         summary        output S from cdfdCycleSummaries
%         diagnostics    diagnostics from cdfdCycleSummaries
%         method         currently 'exact'

    opts = parseCycleAnalysisOptions(varargin{:});

    [C, D, bffInfo] = cdfd_bff(W, ...
        'ToleranceZero', opts.Tol, ...
        'Validate', opts.BFFValidate);

    [Tcycles, cycleStatus] = cdfdExactAllCycles(C, ...
        'Tol', opts.Tol, ...
        'MaxCycles', opts.MaxCycles, ...
        'CheckBalance', opts.CheckBalance, ...
        'CheckReconstruction', opts.CheckReconstruction);

    [summary, diagnostics] = cdfdCycleSummaries(W, C, Tcycles, ...
        'Tol', opts.Tol, ...
        'CheckVolumeIdentity', true, ...
        'CheckReconstruction', opts.CheckReconstruction);

    R = struct();
    R.method = 'exact';
    R.W = W;
    R.C = C;
    R.D = D;
    R.bffInfo = bffInfo;
    R.Tcycles = Tcycles;
    R.cycleStatus = cycleStatus;
    R.summary = summary;
    R.diagnostics = diagnostics;
end


function opts = parseCycleAnalysisOptions(varargin)

    opts = struct();
    opts.Tol = 1e-12;
    opts.MaxCycles = 50000;
    opts.CheckBalance = true;
    opts.CheckReconstruction = true;
    opts.BFFValidate = true;

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

            case 'bffvalidate'
                opts.BFFValidate = value;

            otherwise
                error('Unknown option: %s', name);
        end
    end
end