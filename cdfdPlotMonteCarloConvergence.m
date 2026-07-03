function out = cdfdPlotMonteCarloConvergence(MC, varargin)
% cdfdPlotMonteCarloConvergence
% Plot Monte Carlo convergence diagnostics for CDFD cycle estimation.
%
% PURPOSE
%   Visualise batch convergence from cdfdCycleMonteCarlo.m. The function
%   plots cumulative estimates, or relative errors when exact/reference
%   values are supplied, for:
%
%       T_C      completed circular throughput
%       Lbar_C   mean cycle length
%
%   against cumulative Markov transitions.
%
% INPUT
%   MC   output struct from cdfdCycleMonteCarlo, with:
%          MC.Convergence.Trace
%
% OPTIONS
%   'PlotMode'           default 'estimate'. One of:
%                           'estimate'
%                           'relative_error'
%   'ExactT_C'           default []. Optional exact/reference T_C.
%   'ExactLbar_C'        default []. Optional exact/reference Lbar_C.
%   'ShowBatchEstimates' default true. Overlay batch-only estimates.
%   'FigureVisible'      default true.
%   'ExportPath'         default ''. If non-empty, export the figure.
%   'ExportContentType'  default 'vector'.
%   'Title'              default ''. Figure title.
%   'Parent'             default []. Existing tiledlayout or axes. Usually [].
%
% OUTPUT
%   out   struct with fields:
%           Figure
%           TiledLayout
%           AxesT
%           AxesLbar
%           Lines
%           Trace
%           PlotData
%           ExportPath
%           Options
%
% EXAMPLE
%   out = cdfdPlotMonteCarloConvergence(MC);
%
%   out = cdfdPlotMonteCarloConvergence(MC, ...
%       'ExactT_C', R.summary.T_C, ...
%       'ExactLbar_C', R.summary.Lbar_C, ...
%       'PlotMode', 'relative_error');

    opts = parsePlotOptions(varargin{:});

    validateMonteCarloInput(MC);

    Trace = MC.Convergence.Trace;

    if isempty(Trace) || height(Trace) == 0
        error(['MC.Convergence.Trace is empty. Re-run cdfdCycleMonteCarlo ' ...
               'with BudgetMode=''transitions'' and non-empty BatchTransitions.']);
    end

    required = { ...
        'CumulativeTransitions', ...
        'T_C', ...
        'Lbar_C', ...
        'Batch_T_C', ...
        'Batch_Lbar_C' ...
    };

    for k = 1:numel(required)
        if ~ismember(required{k}, Trace.Properties.VariableNames)
            error('MC.Convergence.Trace is missing required variable: %s.', required{k});
        end
    end

    x = full(Trace.CumulativeTransitions);
    TC = full(Trace.T_C);
    Lbar = full(Trace.Lbar_C);

    BatchTC = full(Trace.Batch_T_C);
    BatchLbar = full(Trace.Batch_Lbar_C);

    plotMode = lower(string(opts.PlotMode));

    switch plotMode
        case "estimate"
            yTC = TC;
            yBatchTC = BatchTC;
            yRefTC = opts.ExactT_C;

            yLbar = Lbar;
            yBatchLbar = BatchLbar;
            yRefLbar = opts.ExactLbar_C;

            yLabelTC = 'Completed circular throughput, $T_C$';
            yLabelLbar = 'Mean cycle length, $\bar{L}_C$';

            %yLabelTC = 'Completed circular throughput, T_C';
            %yLabelLbar = 'Mean cycle length, Lbar_C';

            titleTC = 'Convergence of completed circular throughput';
            titleLbar = 'Convergence of mean cycle length';

            referenceLabel = 'Reference value';

        case "relative_error"
            validateReferenceValue(opts.ExactT_C, 'ExactT_C');
            validateReferenceValue(opts.ExactLbar_C, 'ExactLbar_C');

            yTC = 100 * (TC - opts.ExactT_C) ./ opts.ExactT_C;
            yBatchTC = 100 * (BatchTC - opts.ExactT_C) ./ opts.ExactT_C;
            yRefTC = 0;

            yLbar = 100 * (Lbar - opts.ExactLbar_C) ./ opts.ExactLbar_C;
            yBatchLbar = 100 * (BatchLbar - opts.ExactLbar_C) ./ opts.ExactLbar_C;
            yRefLbar = 0;

            yLabelTC = 'Relative error in $T_C$ (\%)';
            yLabelLbar = 'Relative error in $\bar{L}_C$ (\%)';

            %yLabelTC = 'Relative error in T_C (%)';
            %yLabelLbar = 'Relative error in \bar{L}_C (%)';

            titleTC = 'Relative error in completed circular throughput';
            titleLbar = 'Relative error in mean cycle length';

            referenceLabel = 'Zero error';

        otherwise
            error('PlotMode must be ''estimate'' or ''relative_error''.');
    end

    if isempty(opts.Parent)
        visible = 'on';
        if ~opts.FigureVisible
            visible = 'off';
        end

        fig = figure( ...
            'Color', 'w', ...
            'Name', 'CDFD Monte Carlo convergence', ...
            'Visible', visible);

        fig.Position = [100 100 1150 720];

        tl = tiledlayout(fig, 2, 1, ...
            'TileSpacing', 'compact', ...
            'Padding', 'compact');

        axT = nexttile(tl, 1);
        axL = nexttile(tl, 2);

    else
        if isa(opts.Parent, 'matlab.graphics.layout.TiledChartLayout')
            tl = opts.Parent;
            fig = ancestor(tl, 'figure');
            axT = nexttile(tl);
            axL = nexttile(tl);

        elseif isa(opts.Parent, 'matlab.graphics.axis.Axes')
            axT = opts.Parent;
            fig = ancestor(axT, 'figure');
            tl = [];
            axL = axes(fig); %#ok<LAXES>
            warning(['Parent supplied as a single axes. Created a second axes; ' ...
                     'using a fresh figure is usually cleaner.']);

        else
            error('Parent must be empty, an axes handle, or a tiledlayout handle.');
        end
    end

    lines = struct();

    %% Top panel: T_C

    hold(axT, 'on');

    if opts.ShowBatchEstimates
        lines.BatchT_C = plot(axT, x, yBatchTC, ...
            ':', ...
            'LineWidth', 1.0, ...
            'DisplayName', 'Batch estimate');
    else
        lines.BatchT_C = gobjects(0);
    end

    lines.T_C = plot(axT, x, yTC, ...
        '-o', ...
        'LineWidth', 1.6, ...
        'MarkerSize', 4, ...
        'DisplayName', 'Cumulative estimate');

    if ~isempty(yRefTC)
        lines.ReferenceT_C = yline(axT, yRefTC, ...
            '--', ...
            'LineWidth', 1.2, ...
            'DisplayName', referenceLabel);
    else
        lines.ReferenceT_C = gobjects(0);
    end

    hold(axT, 'off');

    grid(axT, 'on');
    xlabel(axT, 'Cumulative Markov transitions');
    ylabel(axT, yLabelTC, 'Interpreter', 'latex');

    title(axT, titleTC, 'Interpreter', 'none');

    legend(axT, 'Location', 'best');

    cleanAxes(axT);

    if plotMode == "relative_error"
        setSymmetricYLimits(axT);
        setSymmetricYLimits(axL);
    end

    %% Bottom panel: Lbar_C

    hold(axL, 'on');

    if opts.ShowBatchEstimates
        lines.BatchLbar_C = plot(axL, x, yBatchLbar, ...
            ':', ...
            'LineWidth', 1.0, ...
            'DisplayName', 'Batch estimate');
    else
        lines.BatchLbar_C = gobjects(0);
    end

    lines.Lbar_C = plot(axL, x, yLbar, ...
        '-o', ...
        'LineWidth', 1.6, ...
        'MarkerSize', 4, ...
        'DisplayName', 'Cumulative estimate');

    if ~isempty(yRefLbar)
        lines.ReferenceLbar_C = yline(axL, yRefLbar, ...
            '--', ...
            'LineWidth', 1.2, ...
            'DisplayName', referenceLabel);
    else
        lines.ReferenceLbar_C = gobjects(0);
    end

    hold(axL, 'off');

    grid(axL, 'on');
    xlabel(axL, 'Cumulative Markov transitions');
    ylabel(axL, yLabelLbar, 'Interpreter', 'latex');

    title(axL, titleLbar, 'Interpreter', 'none');

    legend(axL, 'Location', 'best');

    cleanAxes(axL);

     if plotMode == "relative_error"
        setSymmetricYLimits(axT);
        setSymmetricYLimits(axL);
    end

    %% Overall title

    if ~isempty(tl)
        if isempty(opts.Title)
            title(tl, 'Monte Carlo convergence diagnostics', ...
                'Interpreter', 'none');
        else
            title(tl, opts.Title, ...
                'Interpreter', 'none');
        end
    elseif ~isempty(opts.Title)
        sgtitle(fig, opts.Title, ...
            'Interpreter', 'none');
    end

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

    PlotData = table( ...
        x, ...
        yTC, ...
        yBatchTC, ...
        yLbar, ...
        yBatchLbar, ...
        'VariableNames', { ...
            'CumulativeTransitions', ...
            'T_C', ...
            'Batch_T_C', ...
            'Lbar_C', ...
            'Batch_Lbar_C' ...
        } ...
    );

    out = struct();
    out.Figure = fig;
    out.TiledLayout = tl;
    out.AxesT = axT;
    out.AxesLbar = axL;
    out.Lines = lines;
    out.Trace = Trace;
    out.PlotData = PlotData;
    out.ExportPath = exportPath;
    out.Options = opts;
end


function validateMonteCarloInput(MC)

    if ~isstruct(MC)
        error('MC must be a struct returned by cdfdCycleMonteCarlo.');
    end

    if ~isfield(MC, 'Convergence')
        error('MC is missing field Convergence.');
    end

    if ~isfield(MC.Convergence, 'Trace')
        error('MC.Convergence is missing field Trace.');
    end

    if ~istable(MC.Convergence.Trace)
        error('MC.Convergence.Trace must be a table.');
    end
end


function validateReferenceValue(x, name)

    if isempty(x)
        error('PlotMode=''relative_error'' requires %s.', name);
    end

    if ~isnumeric(x) || ~isscalar(x) || ~isfinite(x)
        error('%s must be a finite numeric scalar.', name);
    end

    if abs(x) <= eps
        error('%s must be non-zero for relative-error plotting.', name);
    end
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


function opts = parsePlotOptions(varargin)

    opts = struct();
    opts.PlotMode = 'estimate';
    opts.ExactT_C = [];
    opts.ExactLbar_C = [];
    opts.ShowBatchEstimates = true;
    opts.FigureVisible = true;
    opts.ExportPath = '';
    opts.ExportContentType = 'vector';
    opts.Title = '';
    opts.Parent = [];

    if mod(numel(varargin), 2) ~= 0
        error('Options must be name-value pairs.');
    end

    for k = 1:2:numel(varargin)
        name = varargin{k};
        value = varargin{k+1};

        switch lower(name)
            case 'plotmode'
                opts.PlotMode = value;

            case 'exactt_c'
                opts.ExactT_C = value;

            case 'exactlbar_c'
                opts.ExactLbar_C = value;

            case 'showbatchestimates'
                opts.ShowBatchEstimates = value;

            case 'figurevisible'
                opts.FigureVisible = value;

            case 'exportpath'
                opts.ExportPath = value;

            case 'exportcontenttype'
                opts.ExportContentType = value;

            case 'title'
                opts.Title = value;

            case 'parent'
                opts.Parent = value;

            otherwise
                error('Unknown option: %s.', name);
        end
    end
end

function setSymmetricYLimits(ax)
% setSymmetricYLimits
% Use symmetric y-limits around zero for relative-error plots.

    yl = ylim(ax);
    m = max(abs(yl));

    if m == 0 || ~isfinite(m)
        m = 1;
    end

    ylim(ax, [-m, m]);
end