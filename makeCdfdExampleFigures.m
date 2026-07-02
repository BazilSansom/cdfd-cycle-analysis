function results = makeCdfdExampleFigures(outputDir, varargin)
% makeCdfdExampleFigures
% Generate validation/example figures for the CDFD cycle-analysis paper.
%
% PURPOSE
%   Produces three example figures using the current CDFD cycle-analysis
%   toolbox pipeline:
%
%       W -> cdfdCycleAnalysis(W)
%         -> cdfd_bff(W)
%         -> cdfdExactAllCycles(C)
%         -> cdfdCycleSummaries(W,C,Tcycles)
%
%   Figures:
%       fig_bidirected_triangle.pdf
%       fig_subdivision_invariance.pdf
%       fig_signed_overlap.pdf
%
% INPUT
%   outputDir   folder for exported figures. Default: pwd.
%
% OPTIONS
%   'Tol'           default 1e-12
%   'MaxCycles'     default 50000
%   'Export'        default true
%   'ShowFigures'   default true
%   'Verbose'       default true
%   'Formats'       default {'pdf'}
%   'BFFValidate'   default true
%
% OUTPUT
%   results     struct containing networks, cycle tables, summaries,
%               diagnostics, signed summary, figure handles, and exported
%               file paths.
%
% NOTES
%   Requires:
%       cdfd_bff.m
%       cdfdCycleAnalysis.m
%       cdfdExactAllCycles.m
%       cdfdCycleSummaries.m
%
%   The signed example currently computes signed summaries locally from the
%   cycle table. This can later be replaced by cdfdSignedCycleSummaries.m.

    if nargin < 1 || isempty(outputDir)
        outputDir = pwd;
    end

    opts = parseExampleFigureOptions(varargin{:});

    if ~isfolder(outputDir)
        mkdir(outputDir);
    end

    visibility = 'on';
    if ~opts.ShowFigures
        visibility = 'off';
    end

    results = struct();
    results.OutputDir = outputDir;
    results.Options = opts;
    results.ExportedFiles = strings(0, 1);

    %% Example 1: unit bidirected triangle

    Wtri = zeros(3);
    Wtri(1,2) = 1; Wtri(2,1) = 1;
    Wtri(2,3) = 1; Wtri(3,2) = 1;
    Wtri(1,3) = 1; Wtri(3,1) = 1;

    Rtri = cdfdCycleAnalysis(Wtri, ...
        'Tol', opts.Tol, ...
        'MaxCycles', opts.MaxCycles, ...
        'BFFValidate', opts.BFFValidate);

    if opts.Verbose
        fprintf('\nBidirected triangle cycle table:\n');
        disp(Rtri.Tcycles);
        disp(Rtri.cycleStatus);

        fprintf('\nBidirected triangle summary:\n');
        fprintf('V_C   = %.10f\n', Rtri.summary.V_C);
        fprintf('T_C   = %.10f\n', Rtri.summary.T_C);
        fprintf('LbarC = %.10f\n', Rtri.summary.Lbar_C);
    end

    fig1 = makeBidirectedTriangleFigure(Wtri, Rtri, visibility);

    exported = exportExampleFigure(fig1, outputDir, ...
        'fig_bidirected_triangle', opts);

    results.ExportedFiles = [results.ExportedFiles; exported(:)];

    results.BidirectedTriangle = Rtri;
    results.Figures.BidirectedTriangle = fig1;

    %% Example 2: subdivision invariance

    Worig = zeros(3);
    Worig(1,2) = 2;
    Worig(2,3) = 1;
    Worig(3,1) = 1;

    Wsub = zeros(4);
    Wsub(1,4) = 2;
    Wsub(4,2) = 2;
    Wsub(2,3) = 1;
    Wsub(3,1) = 1;

    Rorig = cdfdCycleAnalysis(Worig, ...
        'Tol', opts.Tol, ...
        'MaxCycles', opts.MaxCycles, ...
        'BFFValidate', opts.BFFValidate);

    Rsub = cdfdCycleAnalysis(Wsub, ...
        'Tol', opts.Tol, ...
        'MaxCycles', opts.MaxCycles, ...
        'BFFValidate', opts.BFFValidate);

    if opts.Verbose
        fprintf('\nOriginal subdivision-example cycle table:\n');
        disp(Rorig.Tcycles);
        disp(Rorig.cycleStatus);

        fprintf('\nSubdivided cycle table:\n');
        disp(Rsub.Tcycles);
        disp(Rsub.cycleStatus);

        fprintf('\nOriginal network:\n');
        printCoreSummary(Rorig.summary);

        fprintf('\nSubdivided network:\n');
        printCoreSummary(Rsub.summary);
    end

    fig2 = makeSubdivisionFigure(Worig, Rorig, Wsub, Rsub, visibility);

    exported = exportExampleFigure(fig2, outputDir, ...
        'fig_subdivision_invariance', opts);

    results.ExportedFiles = [results.ExportedFiles; exported(:)];

    results.Subdivision.Original = Rorig;
    results.Subdivision.Subdivided = Rsub;
    results.Figures.Subdivision = fig2;

    %% Example 3: signed overlap and contextual edge roles

    Sigma = zeros(3);

    % Positive signs on all supported edges by default.
    Sigma(Wtri > 0) = 1;

    % Negative signs on the bidirectional pair 1 <-> 2.
    Sigma(1,2) = -1;
    Sigma(2,1) = -1;

    signed = localSignedCycleSummaries(Rtri.Tcycles, Sigma, [1 2]);

    if opts.Verbose
        fprintf('\nSigned summary:\n');
        fprintf('T+        = %.10f\n', signed.Tplus);
        fprintf('T-        = %.10f\n', signed.Tminus);
        fprintf('T_C       = %.10f\n', signed.Ttotal);
        fprintf('V+        = %.10f\n', signed.Vplus);
        fprintf('V-        = %.10f\n', signed.Vminus);
        fprintf('V_C       = %.10f\n', signed.Vtotal);
        fprintf('T+ share  = %.10f\n', signed.ShareTplus);
        fprintf('T- share  = %.10f\n', signed.ShareTminus);
        fprintf('V+ share  = %.10f\n', signed.ShareVplus);
        fprintf('V- share  = %.10f\n', signed.ShareVminus);

        fprintf('\nContextual role of edge 1 -> 2:\n');
        fprintf('c12+      = %.10f\n', signed.TargetEdgeCplus);
        fprintf('c12-      = %.10f\n', signed.TargetEdgeCminus);
        fprintf('c12 total = %.10f\n', signed.TargetEdgeCtotal);
    end

    fig3 = makeSignedOverlapFigure(Wtri, Sigma, signed, visibility);

    exported = exportExampleFigure(fig3, outputDir, ...
        'fig_signed_overlap', opts);

    results.ExportedFiles = [results.ExportedFiles; exported(:)];

    results.Signed.Sigma = Sigma;
    results.Signed.Summary = signed;
    results.Figures.SignedOverlap = fig3;

    if opts.Verbose && opts.Export
        fprintf('\nExported files:\n');
        disp(results.ExportedFiles);
    end
end


%% Figure builders

function fig = makeBidirectedTriangleFigure(Wtri, Rtri, visibility)

    Ttri = Rtri.Tcycles;
    S = Rtri.summary;

    plotTable = Ttri;
    plotTable.Label = makeCycleLabels(plotTable);
    plotTable = sortrows(plotTable, {'Length', 'Label'});

    fig = figure( ...
        'Color', 'w', ...
        'Name', 'Bidirected triangle exact cycle weights', ...
        'Visible', visibility);

    fig.Position = [100 100 1250 650];

    tl = tiledlayout(fig, 1, 2, ...
        'TileSpacing', 'compact', ...
        'Padding', 'compact');

    % Left panel: network.
    ax1 = nexttile(tl, 1);

    G = digraph(Wtri > 0);

    plot(G, ...
        'Layout', 'circle', ...
        'NodeLabel', {'1','2','3'}, ...
        'LineWidth', 1.5, ...
        'ArrowSize', 12);

    title(ax1, 'Bidirected triangle');

    subtitle(ax1, sprintf('$V_C = %.0f,\\ T_C = %.3f,\\ \\bar{L}_C = %.3f$', ...
        S.V_C, S.T_C, S.Lbar_C), ...
        'Interpreter', 'latex');

    axis(ax1, 'off');
    cleanAxes(ax1);

    text(ax1, 0.5, -0.08, 'All directed edge weights = 1', ...
        'Units', 'normalized', ...
        'HorizontalAlignment', 'center', ...
        'Interpreter', 'none');

    % Right panel: exact cycle weights.
    ax2 = nexttile(tl, 2);

    barh(ax2, plotTable.lambda);

    set(ax2, ...
        'YTick', 1:height(plotTable), ...
        'YTickLabel', plotTable.Label, ...
        'YDir', 'reverse');

    ax2.TickLabelInterpreter = 'none';

    xlabel(ax2, '$\lambda_\gamma$', 'Interpreter', 'latex');
    title(ax2, 'Exact probabilistic cycle weights');
    grid(ax2, 'on');
    cleanAxes(ax2);

    xmax = max(plotTable.lambda);
    xlim(ax2, [0, 1.25 * xmax]);

    for r = 1:height(plotTable)
        text(ax2, plotTable.lambda(r), r, ...
            sprintf('  %.3f', plotTable.lambda(r)), ...
            'VerticalAlignment', 'middle');
    end
end


function fig = makeSubdivisionFigure(Worig, Rorig, Wsub, Rsub, visibility)

    Sorig = Rorig.summary;
    Ssub = Rsub.summary;

    fig = figure( ...
        'Color', 'w', ...
        'Name', 'Subdivision invariance example', ...
        'Visible', visibility);

    fig.Position = [100 100 1250 650];

    tl = tiledlayout(fig, 1, 2, ...
        'TileSpacing', 'compact', ...
        'Padding', 'compact');

    % Left panel: original network.
    ax1 = nexttile(tl, 1);

    Gorig = digraph(Worig > 0);

    p = plot(Gorig, ...
        'XData', [0, 1, 0.5], ...
        'YData', [0, 0, 0.85], ...
        'NodeLabel', {'1','2','3'}, ...
        'LineWidth', 1.8, ...
        'ArrowSize', 13);

    p.EdgeLabel = makeEdgeLabelsWC(Worig, Rorig.C, Gorig);

    title(ax1, '(a) Original network');

    subtitle(ax1, sprintf('$\\Gamma_{vol}=%.3f,\\ \\Gamma_{thr}=%.3f$', ...
        Sorig.Gamma_vol, Sorig.Gamma_thr), ...
        'Interpreter', 'latex');

    axis(ax1, 'off');
    cleanAxes(ax1);

    % Right panel: subdivided network.
    ax2 = nexttile(tl, 2);

    Gsub = digraph(Wsub > 0);

    pp = plot(Gsub, ...
        'XData', [0, 1, 0.5, 0.5], ...
        'YData', [0, 0, 0.85, 0], ...
        'NodeLabel', {'1','2','3','4'}, ...
        'LineWidth', 1.8, ...
        'ArrowSize', 13);

    pp.EdgeLabel = makeEdgeLabelsWC(Wsub, Rsub.C, Gsub);

    title(ax2, '(b) Edge 1 -> 2 subdivided');

    subtitle(ax2, sprintf('$\\Gamma_{vol}=%.3f,\\ \\Gamma_{thr}=%.3f$', ...
        Ssub.Gamma_vol, Ssub.Gamma_thr), ...
        'Interpreter', 'latex');

    axis(ax2, 'off');
    cleanAxes(ax2);

    annotation(fig, 'textbox', [0.32, 0.01, 0.36, 0.05], ...
        'String', 'Edge labels show (w,c): total flow and circular flow.', ...
        'EdgeColor', 'none', ...
        'HorizontalAlignment', 'center', ...
        'Interpreter', 'none');
end


function fig = makeSignedOverlapFigure(Wtri, Sigma, signed, visibility)

    plotTable = signed.CycleTable;
    plotTable.Label = makeCycleLabels(plotTable);

    plotTable.PolarityLabel = strings(height(plotTable), 1);
    plotTable.PolarityLabel(plotTable.Polarity == 1) = "positive";
    plotTable.PolarityLabel(plotTable.Polarity == -1) = "negative";

    % Positive cycles first.
    plotTable.SortPolarity = -plotTable.Polarity;
    plotTable = sortrows(plotTable, {'SortPolarity', 'Length', 'Label'});

    fig = figure( ...
        'Color', 'w', ...
        'Name', 'Signed overlap and contextual edge roles', ...
        'Visible', visibility);

    fig.Position = [100 100 1250 650];

    tl = tiledlayout(fig, 1, 2, ...
        'TileSpacing', 'compact', ...
        'Padding', 'compact');

    % Left panel: signed magnitude network.
    ax1 = nexttile(tl, 1);

    G = digraph(Wtri > 0);

    xNode = [1.0, 0.0, 0.0];
    yNode = [0.0, 0.8, -0.8];

    p = plot(G, ...
        'XData', xNode, ...
        'YData', yNode, ...
        'NodeLabel', {'1','2','3'}, ...
        'LineWidth', 1.4, ...
        'ArrowSize', 12);

    p.EdgeColor = [0.45 0.70 0.90];
    p.EdgeLabel = repmat("", numedges(G), 1);

    try
        p.NodeFontSize = 13;
    catch
    end

    highlight(p, 1, 2, ...
        'EdgeColor', [0.20 0.20 0.20], ...
        'LineWidth', 4.5, ...
        'ArrowSize', 18);

    addManualSignLabels(ax1, Sigma);

    title(ax1, '(a) Signed magnitude network');

    subtitle(ax1, sprintf('$c_{12}^{+}=%.3f,\\ c_{12}^{-}=%.3f,\\ c_{12}=%.3f$', ...
        signed.TargetEdgeCplus, ...
        signed.TargetEdgeCminus, ...
        signed.TargetEdgeCtotal), ...
        'Interpreter', 'latex');

    axis(ax1, 'off');
    xlim(ax1, [-0.20, 1.15]);
    ylim(ax1, [-0.95, 0.95]);
    cleanAxes(ax1);

    text(ax1, 0.5, -0.08, ...
        'Bold edge is 1 -> 2; edge labels show causal signs.', ...
        'Units', 'normalized', ...
        'HorizontalAlignment', 'center', ...
        'Interpreter', 'none');

    % Right panel: cycle weights by polarity.
    ax2 = nexttile(tl, 2);

    b = barh(ax2, plotTable.lambda);
    b.FaceColor = 'flat';

    for r = 1:height(plotTable)
        if plotTable.Polarity(r) == 1
            b.CData(r, :) = [0.25, 0.55, 0.25];
        else
            b.CData(r, :) = [0.75, 0.25, 0.25];
        end
    end

    set(ax2, ...
        'YTick', 1:height(plotTable), ...
        'YTickLabel', plotTable.Label, ...
        'YDir', 'reverse');

    ax2.TickLabelInterpreter = 'none';

    xlabel(ax2, '$\lambda_\gamma$', 'Interpreter', 'latex');
    title(ax2, '(b) Cycle weights by polarity');
    grid(ax2, 'on');
    cleanAxes(ax2);

    xmax = max(plotTable.lambda);
    xlim(ax2, [0, 1.35 * xmax]);

    for r = 1:height(plotTable)
        if plotTable.Polarity(r) == 1
            polLabel = '+';
        else
            polLabel = '-';
        end

        text(ax2, plotTable.lambda(r), r, ...
            sprintf('  %.3f (%s)', plotTable.lambda(r), polLabel), ...
            'VerticalAlignment', 'middle');
    end

    subtitle(ax2, sprintf('$T^{+}=%.3f,\\ T^{-}=%.3f,\\ V^{+}=%.3f,\\ V^{-}=%.3f$', ...
        signed.Tplus, signed.Tminus, signed.Vplus, signed.Vminus), ...
        'Interpreter', 'latex');
end


%% Signed local summary

function signed = localSignedCycleSummaries(Tcycles, Sigma, targetEdge)
% localSignedCycleSummaries
% Compute signed feedback summaries from a cycle table and sign matrix.
%
% This is a local helper for the example figure. It can later be promoted
% into a separate cdfdSignedCycleSummaries.m utility.

    Tsigned = Tcycles;

    polarity = zeros(height(Tsigned), 1);
    containsTarget = false(height(Tsigned), 1);

    for r = 1:height(Tsigned)
        cyc = Tsigned.Cycle{r};
        closedCyc = [cyc(:).' cyc(1)];

        sig = 1;

        for k = 1:numel(cyc)
            i = closedCyc(k);
            j = closedCyc(k+1);

            sig = sig * Sigma(i,j);

            if i == targetEdge(1) && j == targetEdge(2)
                containsTarget(r) = true;
            end
        end

        polarity(r) = sig;
    end

    Tsigned.Polarity = polarity;
    Tsigned.ContainsTargetEdge = containsTarget;

    positive = Tsigned.Polarity == 1;
    negative = Tsigned.Polarity == -1;

    signed = struct();
    signed.CycleTable = Tsigned;

    signed.Tplus = sum(Tsigned.lambda(positive));
    signed.Tminus = sum(Tsigned.lambda(negative));
    signed.Ttotal = signed.Tplus + signed.Tminus;

    signed.Vplus = sum(Tsigned.lambda(positive) .* Tsigned.Length(positive));
    signed.Vminus = sum(Tsigned.lambda(negative) .* Tsigned.Length(negative));
    signed.Vtotal = signed.Vplus + signed.Vminus;

    signed.ShareTplus = signed.Tplus / signed.Ttotal;
    signed.ShareTminus = signed.Tminus / signed.Ttotal;

    signed.ShareVplus = signed.Vplus / signed.Vtotal;
    signed.ShareVminus = signed.Vminus / signed.Vtotal;

    signed.TargetEdge = targetEdge;
    signed.TargetEdgeCplus = sum(Tsigned.lambda(positive & containsTarget));
    signed.TargetEdgeCminus = sum(Tsigned.lambda(negative & containsTarget));
    signed.TargetEdgeCtotal = signed.TargetEdgeCplus + signed.TargetEdgeCminus;
end


%% Export and display helpers

function exported = exportExampleFigure(fig, outputDir, baseName, opts)

    exported = strings(0, 1);

    if ~opts.Export
        return;
    end

    formats = string(opts.Formats);

    for k = 1:numel(formats)
        fmt = lower(formats(k));
        filename = fullfile(outputDir, baseName + "." + fmt);

        switch fmt
            case "pdf"
                exportgraphics(fig, filename, 'ContentType', 'vector');

            case {"png", "jpg", "jpeg", "tif", "tiff"}
                exportgraphics(fig, filename, 'Resolution', 300);

            otherwise
                error('Unsupported export format: %s', fmt);
        end

        exported(end+1, 1) = filename; %#ok<AGROW>
        fprintf('Exported %s\n', filename);
    end
end


function printCoreSummary(S)

    fprintf('V_C       = %.10f\n', S.V_C);
    fprintf('V_W       = %.10f\n', S.V_W);
    fprintf('T_C       = %.10f\n', S.T_C);
    fprintf('D_net     = %.10f\n', S.D_net);
    fprintf('Gamma_vol = %.10f\n', S.Gamma_vol);
    fprintf('Gamma_thr = %.10f\n', S.Gamma_thr);
end


function labels = makeCycleLabels(T)
% makeCycleLabels
% Convert cycle node lists to readable labels.
% Self-loops are shown as "i -> i".
% Two-cycles are shown as "i <-> j".
% Longer cycles show closure.

    labels = strings(height(T), 1);

    for r = 1:height(T)
        cyc = T.Cycle{r};

        if numel(cyc) == 1
            labels(r) = sprintf('%d -> %d', cyc(1), cyc(1));

        elseif numel(cyc) == 2
            labels(r) = sprintf('%d <-> %d', cyc(1), cyc(2));

        else
            closedCyc = [cyc(:).' cyc(1)];
            labels(r) = strjoin(string(closedCyc), " -> ");
        end
    end
end


function labels = makeEdgeLabelsWC(W, C, G)
% makeEdgeLabelsWC
% Edge labels of the form "(w,c)", where w is total flow and c is circular flow.

    W = full(W);
    C = full(C);

    E = G.Edges.EndNodes;
    labels = strings(size(E, 1), 1);

    for r = 1:size(E, 1)
        i = E(r, 1);
        j = E(r, 2);
        labels(r) = sprintf('(%g,%g)', W(i,j), C(i,j));
    end
end


function addManualSignLabels(ax, Sigma)
% addManualSignLabels
% Manual sign labels for the specific signed bidirected triangle layout.
%
% Node coordinates are:
%   1: ( 1.0,  0.0)
%   2: ( 0.0,  0.8)
%   3: ( 0.0, -0.8)

    %#ok<INUSD>
    % These locations are chosen visually for the fixed three-node example.

    text(ax, 0.56, 0.50, '-', ...
        'FontSize', 16, ...
        'FontWeight', 'bold', ...
        'HorizontalAlignment', 'center', ...
        'Interpreter', 'none');

    text(ax, 0.42, 0.62, '-', ...
        'FontSize', 16, ...
        'FontWeight', 'bold', ...
        'HorizontalAlignment', 'center', ...
        'Interpreter', 'none');

    text(ax, -0.055, 0.10, '+', ...
        'FontSize', 16, ...
        'FontWeight', 'bold', ...
        'HorizontalAlignment', 'center', ...
        'Interpreter', 'none');

    text(ax, 0.105, -0.10, '+', ...
        'FontSize', 16, ...
        'FontWeight', 'bold', ...
        'HorizontalAlignment', 'center', ...
        'Interpreter', 'none');

    text(ax, 0.52, -0.42, '+', ...
        'FontSize', 16, ...
        'FontWeight', 'bold', ...
        'HorizontalAlignment', 'center', ...
        'Interpreter', 'none');

    text(ax, 0.38, -0.58, '+', ...
        'FontSize', 16, ...
        'FontWeight', 'bold', ...
        'HorizontalAlignment', 'center', ...
        'Interpreter', 'none');
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


function opts = parseExampleFigureOptions(varargin)

    opts = struct();
    opts.Tol = 1e-12;
    opts.MaxCycles = 50000;
    opts.Export = true;
    opts.ShowFigures = true;
    opts.Verbose = true;
    opts.Formats = {'pdf'};
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

            case 'export'
                opts.Export = value;

            case 'showfigures'
                opts.ShowFigures = value;

            case 'verbose'
                opts.Verbose = value;

            case 'formats'
                opts.Formats = value;

            case 'bffvalidate'
                opts.BFFValidate = value;

            otherwise
                error('Unknown option: %s', name);
        end
    end
end