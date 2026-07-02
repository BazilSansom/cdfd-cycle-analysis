function out = cdfdPlotCycleCoverage(Tcycles, varargin)
% cdfdPlotCycleCoverage
% Plot cumulative cycle-coverage curves from probabilistic cycle coefficients.
%
% PURPOSE
%   Visualise how much circular throughput and circular edge-volume are
%   accounted for by the top-ranked cycles.
%
%   The function ranks cycles using cdfdTopCycles and plots:
%
%       cumulative throughput coverage
%       cumulative volume coverage
%
%   where
%
%       I_thr(gamma) = lambda_gamma
%       I_vol(gamma) = lambda_gamma * |gamma|.
%
% INPUT
%   Tcycles   cycle table, usually returned by cdfdExactAllCycles(C), with
%             at least columns:
%               Cycle, Length, lambda
%
% OPTIONS
%   'TopK'              default Inf. Number of ranked cycles to plot.
%   'RankBy'            default 'throughput'. One of:
%                         'throughput'
%                         'volume'
%                         'length'
%   'Tol'               default 1e-12
%   'ShowThroughput'    default true
%   'ShowVolume'        default true
%   'Percent'           default true. If true, plot coverage in percent.
%   'Title'             default ''. Figure title. If empty, a default title
%                       is generated.
%   'FigureVisible'     default true
%   'Parent'            default []. Existing axes to plot into.
%   'ExportPath'        default ''. If non-empty, export the figure.
%   'ExportContentType' default 'vector'
%   'WeightVariable'    default 'lambda'
%   'LengthVariable'    default 'Length'
%   'CycleVariable'     default 'Cycle'
%
% OUTPUT
%   out   struct with fields:
%           Figure
%           Axes
%           TopCycles
%           Lines
%           Summary
%           ExportPath
%           Options
%
% EXAMPLE
%   R = cdfdCycleAnalysis(W);
%   out = cdfdPlotCycleCoverage(R.Tcycles, ...
%       'TopK', 20, ...
%       'RankBy', 'throughput');

    opts = parseCoveragePlotOptions(varargin{:});
    tol = opts.Tol;

    if ~opts.ShowThroughput && ~opts.ShowVolume
        error('At least one of ShowThroughput or ShowVolume must be true.');
    end

    Top = cdfdTopCycles(Tcycles, ...
        'TopK', opts.TopK, ...
        'RankBy', opts.RankBy, ...
        'Tol', opts.Tol, ...
        'WeightVariable', opts.WeightVariable, ...
        'LengthVariable', opts.LengthVariable, ...
        'CycleVariable', opts.CycleVariable);

    if isempty(Top)
        error('Tcycles is empty; no coverage plot can be produced.');
    end

    x = [0; Top.Rank];

    yThroughput = [0; Top.CumulativeThroughputShare];
    yVolume = [0; Top.CumulativeVolumeShare];

    yLabel = 'Cumulative coverage';

    if opts.Percent
        yThroughput = 100 * yThroughput;
        yVolume = 100 * yVolume;
        yLabel = 'Cumulative coverage (%)';
    end

    if isempty(opts.Parent)
        visible = 'on';
        if ~opts.FigureVisible
            visible = 'off';
        end

        fig = figure( ...
            'Color', 'w', ...
            'Name', 'CDFD cycle coverage', ...
            'Visible', visible);

        ax = axes(fig);
    else
        ax = opts.Parent;
        fig = ancestor(ax, 'figure');
    end

    holdState = ishold(ax);
    hold(ax, 'on');

    lineHandles = struct();

    if opts.ShowThroughput
        lineHandles.Throughput = plot(ax, x, yThroughput, ...
            '-o', ...
            'LineWidth', 1.5, ...
            'MarkerSize', 5, ...
            'DisplayName', 'Throughput coverage');
    else
        lineHandles.Throughput = gobjects(0);
    end

    if opts.ShowVolume
        lineHandles.Volume = plot(ax, x, yVolume, ...
            '-s', ...
            'LineWidth', 1.5, ...
            'MarkerSize', 5, ...
            'DisplayName', 'Volume coverage');
    else
        lineHandles.Volume = gobjects(0);
    end

    if ~holdState
        hold(ax, 'off');
    end

    grid(ax, 'on');

    xlabel(ax, 'Number of top cycles, k');
    ylabel(ax, yLabel);

    if opts.Percent
        ylim(ax, [0 105]);
    else
        ylim(ax, [0 1.05]);
    end

    xlim(ax, [0 max(1, max(x))]);
    xticks(ax, 0:max(x));

    if isempty(opts.Title)
        titleText = sprintf('Cycle coverage ranked by %s', lower(string(opts.RankBy)));
    else
        titleText = opts.Title;
    end

    title(ax, titleText, 'Interpreter', 'none');

    legend(ax, 'Location', 'southeast');

    cleanAxes(ax);

    summary = makeCoverageSummary(Tcycles, Top, opts, tol);

    exportPath = "";

    if strlength(string(opts.ExportPath)) > 0
        exportPath = string(opts.ExportPath);

        [exportDir, ~, ~] = fileparts(exportPath);
        if strlength(exportDir) > 0 && ~isfolder(exportDir)
            mkdir(exportDir);
        end

        exportgraphics(fig, exportPath, ...
            'ContentType', opts.ExportContentType);

        fprintf('Exported %s\n', exportPath);
    end

    out = struct();
    out.Figure = fig;
    out.Axes = ax;
    out.TopCycles = Top;
    out.Lines = lineHandles;
    out.Summary = summary;
    out.ExportPath = exportPath;
    out.Options = opts;
end


function summary = makeCoverageSummary(Tcycles, Top, opts, tol)
% makeCoverageSummary
% Compute compact scalar summary information for the plot.

    lambda = full(Tcycles.(opts.WeightVariable));
    lengths = full(Tcycles.(opts.LengthVariable));

    lambda = lambda(:);
    lengths = lengths(:);

    lambda(abs(lambda) <= tol) = 0;
    lengths = round(lengths);

    volumeContribution = lambda .* lengths;

    T_C = full(sum(lambda));
    V_C = full(sum(volumeContribution));

    if T_C > tol
        p = lambda / T_C;
        HHI = full(sum(p.^2));
        N_eff = 1 / HHI;
    else
        HHI = NaN;
        N_eff = NaN;
    end

    summary = struct();
    summary.NumCyclesTotal = height(Tcycles);
    summary.NumCyclesPlotted = height(Top);
    summary.RankBy = string(opts.RankBy);
    summary.T_C = T_C;
    summary.V_C = V_C;
    summary.HHI = HHI;
    summary.N_eff = N_eff;
    summary.FinalThroughputCoverage = Top.CumulativeThroughputShare(end);
    summary.FinalVolumeCoverage = Top.CumulativeVolumeShare(end);
end


function cleanAxes(ax)
% cleanAxes
% Remove interactive toolbar artefacts before export where supported.

    try
        disableDefaultInteractivity(ax);
    catch
    end

    try
        ax.Toolbar.Visible = 'off';
    catch
    end
end


function opts = parseCoveragePlotOptions(varargin)

    opts = struct();
    opts.TopK = Inf;
    opts.RankBy = 'throughput';
    opts.Tol = 1e-12;

    opts.ShowThroughput = true;
    opts.ShowVolume = true;
    opts.Percent = true;

    opts.Title = '';
    opts.FigureVisible = true;
    opts.Parent = [];

    opts.ExportPath = '';
    opts.ExportContentType = 'vector';

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

            case 'tol'
                opts.Tol = value;

            case 'showthroughput'
                opts.ShowThroughput = value;

            case 'showvolume'
                opts.ShowVolume = value;

            case 'percent'
                opts.Percent = value;

            case 'title'
                opts.Title = value;

            case 'figurevisible'
                opts.FigureVisible = value;

            case 'parent'
                opts.Parent = value;

            case 'exportpath'
                opts.ExportPath = value;

            case 'exportcontenttype'
                opts.ExportContentType = value;

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