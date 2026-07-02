function results = runBasicValidation(varargin)
% runBasicValidation
% Basic validation tests for the CDFD cycle-analysis toolbox.
%
% USAGE
%   results = runBasicValidation();
%   results = runBasicValidation('Verbose', true);
%
% PURPOSE
%   Runs small deterministic examples with known analytical outputs:
%
%     1. Self-loop
%     2. Directed triangle
%     3. Triangle plus directional excess
%     4. Subdivision invariance
%     5. Bidirected triangle exact cycle weights
%     6. Exact nominated cycle weights
%     7. Signed feedback summaries
%     8. Edge co-participation kernel
%
%   The function uses assert-style checks. If it completes without error,
%   the basic toolbox functionality is consistent with the paper examples.

    opts = parseValidationOptions(varargin{:});
    tol = opts.Tol;

    addRepoRootToPath();

    results = struct();
    results.Passed = false;
    results.Tests = table();

    testNames = strings(0, 1);
    testPassed = false(0, 1);
    testError = strings(0, 1);

    if opts.Verbose
        fprintf('\nRunning CDFD cycle-analysis basic validation...\n');
    end

    %% 1. Self-loop

    testName = "Self-loop";

    try
        Wself = 5;

        Rself = cdfdCycleAnalysis(Wself, 'Tol', tol);

        assertClose(Rself.summary.V_W, 5, tol, 'Self-loop V_W');
        assertClose(Rself.summary.V_C, 5, tol, 'Self-loop V_C');
        assertClose(Rself.summary.T_C, 5, tol, 'Self-loop T_C');
        assertClose(Rself.summary.Lbar_C, 1, tol, 'Self-loop Lbar_C');
        assertClose(Rself.summary.Gamma_vol, 1, tol, 'Self-loop Gamma_vol');
        assertClose(Rself.summary.Gamma_thr, 1, tol, 'Self-loop Gamma_thr');

        assert(height(Rself.Tcycles) == 1, 'Self-loop should have one cycle.');
        assert(Rself.Tcycles.Length(1) == 1, 'Self-loop cycle should have length one.');
        assertClose(Rself.Tcycles.lambda(1), 5, tol, 'Self-loop lambda');

        outSelf = cdfdExactCycleWeight(Rself.C, [1]);
        assertClose(outSelf.lambda, 5, tol, 'Nominated self-loop lambda');
        assertClose(outSelf.q, 1, tol, 'Nominated self-loop q');

        [testNames, testPassed, testError] = recordPass(testNames, testPassed, testError, testName);

    catch ME
        [testNames, testPassed, testError] = recordFail(testNames, testPassed, testError, testName, ME);
        if opts.StopOnFailure
            rethrow(ME);
        end
    end

    %% 2. Directed triangle

    testName = "Directed triangle";

    try
        Wcycle = zeros(3);
        Wcycle(1,2) = 1;
        Wcycle(2,3) = 1;
        Wcycle(3,1) = 1;

        Rcycle = cdfdCycleAnalysis(Wcycle, 'Tol', tol);

        assertClose(Rcycle.summary.V_W, 3, tol, 'Directed triangle V_W');
        assertClose(Rcycle.summary.V_C, 3, tol, 'Directed triangle V_C');
        assertClose(Rcycle.summary.T_C, 1, tol, 'Directed triangle T_C');
        assertClose(Rcycle.summary.Lbar_C, 3, tol, 'Directed triangle Lbar_C');
        assertClose(Rcycle.summary.Gamma_vol, 1, tol, 'Directed triangle Gamma_vol');
        assertClose(Rcycle.summary.Gamma_thr, 1, tol, 'Directed triangle Gamma_thr');

        assert(height(Rcycle.Tcycles) == 1, 'Directed triangle should have one cycle.');
        assertClose(Rcycle.Tcycles.lambda(1), 1, tol, 'Directed triangle lambda');

        [testNames, testPassed, testError] = recordPass(testNames, testPassed, testError, testName);

    catch ME
        [testNames, testPassed, testError] = recordFail(testNames, testPassed, testError, testName, ME);
        if opts.StopOnFailure
            rethrow(ME);
        end
    end

    %% 3. Triangle plus directional excess

    testName = "Triangle plus directional excess";

    try
        Wexcess = zeros(3);
        Wexcess(1,2) = 2;
        Wexcess(2,3) = 1;
        Wexcess(3,1) = 1;

        Rexcess = cdfdCycleAnalysis(Wexcess, 'Tol', tol);

        assertClose(Rexcess.summary.V_W, 4, tol, 'Triangle plus direction V_W');
        assertClose(Rexcess.summary.V_C, 3, tol, 'Triangle plus direction V_C');
        assertClose(Rexcess.summary.T_C, 1, tol, 'Triangle plus direction T_C');
        assertClose(Rexcess.summary.D_net, 1, tol, 'Triangle plus direction D_net');
        assertClose(Rexcess.summary.Gamma_vol, 0.75, tol, 'Triangle plus direction Gamma_vol');
        assertClose(Rexcess.summary.Gamma_thr, 0.5, tol, 'Triangle plus direction Gamma_thr');

        [testNames, testPassed, testError] = recordPass(testNames, testPassed, testError, testName);

    catch ME
        [testNames, testPassed, testError] = recordFail(testNames, testPassed, testError, testName, ME);
        if opts.StopOnFailure
            rethrow(ME);
        end
    end

    %% 4. Subdivision invariance

    testName = "Subdivision invariance";

    try
        Worig = zeros(3);
        Worig(1,2) = 2;
        Worig(2,3) = 1;
        Worig(3,1) = 1;

        Wsub = zeros(4);
        Wsub(1,4) = 2;
        Wsub(4,2) = 2;
        Wsub(2,3) = 1;
        Wsub(3,1) = 1;

        Rorig = cdfdCycleAnalysis(Worig, 'Tol', tol);
        Rsub = cdfdCycleAnalysis(Wsub, 'Tol', tol);

        assertClose(Rorig.summary.T_C, 1, tol, 'Original T_C');
        assertClose(Rsub.summary.T_C, 1, tol, 'Subdivided T_C');

        assertClose(Rorig.summary.D_net, 1, tol, 'Original D_net');
        assertClose(Rsub.summary.D_net, 1, tol, 'Subdivided D_net');

        assertClose(Rorig.summary.Gamma_thr, 0.5, tol, 'Original Gamma_thr');
        assertClose(Rsub.summary.Gamma_thr, 0.5, tol, 'Subdivided Gamma_thr');

        assertClose(Rorig.summary.Gamma_vol, 0.75, tol, 'Original Gamma_vol');
        assertClose(Rsub.summary.Gamma_vol, 2/3, tol, 'Subdivided Gamma_vol');

        assertClose(Rorig.summary.Gamma_thr, Rsub.summary.Gamma_thr, tol, ...
            'Throughput circularity should be subdivision invariant');

        assert(abs(Rorig.summary.Gamma_vol - Rsub.summary.Gamma_vol) > 1e-6, ...
            'Edge-volume circularity should change in this subdivision example.');

        [testNames, testPassed, testError] = recordPass(testNames, testPassed, testError, testName);

    catch ME
        [testNames, testPassed, testError] = recordFail(testNames, testPassed, testError, testName, ME);
        if opts.StopOnFailure
            rethrow(ME);
        end
    end

    %% 5. Bidirected triangle exact cycle weights

    testName = "Bidirected triangle";

    try
        Wtri = zeros(3);
        Wtri(1,2) = 1; Wtri(2,1) = 1;
        Wtri(2,3) = 1; Wtri(3,2) = 1;
        Wtri(1,3) = 1; Wtri(3,1) = 1;

        Rtri = cdfdCycleAnalysis(Wtri, 'Tol', tol);

        assertClose(Rtri.summary.V_W, 6, tol, 'Bidirected triangle V_W');
        assertClose(Rtri.summary.V_C, 6, tol, 'Bidirected triangle V_C');
        assertClose(Rtri.summary.T_C, 8/3, tol, 'Bidirected triangle T_C');
        assertClose(Rtri.summary.Lbar_C, 9/4, tol, 'Bidirected triangle Lbar_C');
        assertClose(Rtri.summary.Gamma_vol, 1, tol, 'Bidirected triangle Gamma_vol');
        assertClose(Rtri.summary.Gamma_thr, 1, tol, 'Bidirected triangle Gamma_thr');

        assert(height(Rtri.Tcycles) == 5, 'Bidirected triangle should have five simple cycles.');

        lambda2 = sort(Rtri.Tcycles.lambda(Rtri.Tcycles.Length == 2));
        lambda3 = sort(Rtri.Tcycles.lambda(Rtri.Tcycles.Length == 3));

        assert(numel(lambda2) == 3, 'Bidirected triangle should have three 2-cycles.');
        assert(numel(lambda3) == 2, 'Bidirected triangle should have two 3-cycles.');

        assertClose(lambda2, repmat(2/3, 3, 1), tol, 'Bidirected triangle 2-cycle lambdas');
        assertClose(lambda3, repmat(1/3, 2, 1), tol, 'Bidirected triangle 3-cycle lambdas');

        results.Rtri = Rtri;

        [testNames, testPassed, testError] = recordPass(testNames, testPassed, testError, testName);

    catch ME
        [testNames, testPassed, testError] = recordFail(testNames, testPassed, testError, testName, ME);
        if opts.StopOnFailure
            rethrow(ME);
        end
    end

    %% 6. Exact nominated cycle weights

    testName = "Exact nominated cycle weights";

    try
        Wtri = zeros(3);
        Wtri(1,2) = 1; Wtri(2,1) = 1;
        Wtri(2,3) = 1; Wtri(3,2) = 1;
        Wtri(1,3) = 1; Wtri(3,1) = 1;

        Rtri = cdfdCycleAnalysis(Wtri, 'Tol', tol);

        out12 = cdfdExactCycleWeight(Rtri.C, [1 2], 'Tol', tol);
        out123 = cdfdExactCycleWeight(Rtri.C, [1 2 3], 'Tol', tol);

        assertClose(out12.lambda, 2/3, tol, 'Nominated 2-cycle lambda');
        assertClose(out12.q, 1/9, tol, 'Nominated 2-cycle q');

        assertClose(out123.lambda, 1/3, tol, 'Nominated 3-cycle lambda');
        assertClose(out123.q, 1/18, tol, 'Nominated 3-cycle q');

        [testNames, testPassed, testError] = recordPass(testNames, testPassed, testError, testName);

    catch ME
        [testNames, testPassed, testError] = recordFail(testNames, testPassed, testError, testName, ME);
        if opts.StopOnFailure
            rethrow(ME);
        end
    end

    %% 7. Signed feedback summaries

    testName = "Signed feedback summaries";

    try
        Wtri = zeros(3);
        Wtri(1,2) = 1; Wtri(2,1) = 1;
        Wtri(2,3) = 1; Wtri(3,2) = 1;
        Wtri(1,3) = 1; Wtri(3,1) = 1;

        Rtri = cdfdCycleAnalysis(Wtri, 'Tol', tol);

        Sigma = zeros(3);
        Sigma(Wtri > 0) = 1;
        Sigma(1,2) = -1;
        Sigma(2,1) = -1;

        signed = cdfdSignedCycleSummaries(Rtri.Tcycles, Sigma, ...
            'TargetEdges', [1 2], ...
            'Tol', tol);

        assertClose(signed.Tplus, 2, tol, 'Signed Tplus');
        assertClose(signed.Tminus, 2/3, tol, 'Signed Tminus');
        assertClose(signed.Ttotal, 8/3, tol, 'Signed Ttotal');

        assertClose(signed.Vplus, 4, tol, 'Signed Vplus');
        assertClose(signed.Vminus, 2, tol, 'Signed Vminus');
        assertClose(signed.Vtotal, 6, tol, 'Signed Vtotal');

        assertClose(signed.ShareTplus, 0.75, tol, 'Signed ShareTplus');
        assertClose(signed.ShareTminus, 0.25, tol, 'Signed ShareTminus');

        assertClose(signed.TargetEdgeCplus, 2/3, tol, 'Signed target edge Cplus');
        assertClose(signed.TargetEdgeCminus, 1/3, tol, 'Signed target edge Cminus');
        assertClose(signed.TargetEdgeCtotal, 1, tol, 'Signed target edge Ctotal');

        results.SignedTriangle = signed;

        [testNames, testPassed, testError] = recordPass(testNames, testPassed, testError, testName);

    catch ME
        [testNames, testPassed, testError] = recordFail(testNames, testPassed, testError, testName, ME);
        if opts.StopOnFailure
            rethrow(ME);
        end
    end

    %% 8. Edge co-participation kernel

    testName = "Co-participation kernel";

    try
        Wtri = zeros(3);
        Wtri(1,2) = 1; Wtri(2,1) = 1;
        Wtri(2,3) = 1; Wtri(3,2) = 1;
        Wtri(1,3) = 1; Wtri(3,1) = 1;

        Rtri = cdfdCycleAnalysis(Wtri, 'Tol', tol);

        Kout = cdfdCoParticipationKernel(Rtri.Tcycles, 3, ...
            'MakePairTable', true, ...
            'Tol', tol);

        diagError = max(abs(diag(Kout.K) - Kout.EdgeTable.Cedge));
        reconError = norm(Kout.CedgeMatrix - full(Rtri.C), 'fro');

        assertClose(diagError, 0, tol, 'diag(K) equals Cedge');
        assertClose(reconError, 0, tol, 'CedgeMatrix reconstructs C');

        assert(issymmetric(Kout.K), 'K should be symmetric.');
        assert(issymmetric(Kout.Khat), 'Khat should be symmetric.');

        assert(all(diag(Kout.Khat) <= 1 + 10*tol), 'diag(Khat) should be at most one.');
        assert(all(diag(Kout.Khat) >= 1 - 10*tol), 'diag(Khat) should be one for positive edges.');

        results.KoutTriangle = Kout;

        [testNames, testPassed, testError] = recordPass(testNames, testPassed, testError, testName);

    catch ME
        [testNames, testPassed, testError] = recordFail(testNames, testPassed, testError, testName, ME);
        if opts.StopOnFailure
            rethrow(ME);
        end
    end

    %% Finish

    results.Tests = table(testNames, testPassed, testError, ...
        'VariableNames', {'Test', 'Passed', 'ErrorMessage'});

    results.Passed = all(testPassed);
    results.NumTests = numel(testPassed);
    results.NumPassed = sum(testPassed);
    results.NumFailed = sum(~testPassed);

    if opts.Verbose
        fprintf('\nValidation results:\n');
        disp(results.Tests);

        if results.Passed
            fprintf('All %d validation tests passed.\n', results.NumTests);
        else
            fprintf('%d of %d validation tests failed.\n', ...
                results.NumFailed, results.NumTests);
        end
    end

    if ~results.Passed && opts.ErrorOnFailure
        error('Basic validation failed. See results.Tests for details.');
    end
end


%% Helper functions

function addRepoRootToPath()
% addRepoRootToPath
% Allows this file to live either in the repository root or in tests/.

    thisFile = mfilename('fullpath');
    thisDir = fileparts(thisFile);

    if exist(fullfile(thisDir, 'cdfd_bff.m'), 'file')
        repoRoot = thisDir;
    elseif exist(fullfile(fileparts(thisDir), 'cdfd_bff.m'), 'file')
        repoRoot = fileparts(thisDir);
    else
        repoRoot = pwd;
    end

    addpath(repoRoot);
end


function assertClose(actual, expected, tol, label)

    if nargin < 4
        label = 'value';
    end

    actual = full(actual);
    expected = full(expected);

    err = max(abs(actual(:) - expected(:)));

    scale = max(1, max(abs(expected(:))));

    if err > 100 * tol * scale
        error('%s mismatch. Max abs error %.3g exceeds tolerance %.3g.', ...
            label, err, 100 * tol * scale);
    end
end


function [names, passed, errors] = recordPass(names, passed, errors, testName)

    names(end+1, 1) = testName;
    passed(end+1, 1) = true;
    errors(end+1, 1) = "";
end


function [names, passed, errors] = recordFail(names, passed, errors, testName, ME)

    names(end+1, 1) = testName;
    passed(end+1, 1) = false;
    errors(end+1, 1) = string(ME.message);
end


function opts = parseValidationOptions(varargin)

    opts = struct();
    opts.Tol = 1e-12;
    opts.Verbose = true;
    opts.StopOnFailure = true;
    opts.ErrorOnFailure = true;

    if mod(numel(varargin), 2) ~= 0
        error('Options must be name-value pairs.');
    end

    for k = 1:2:numel(varargin)
        name = varargin{k};
        value = varargin{k+1};

        switch lower(name)
            case 'tol'
                opts.Tol = value;

            case 'verbose'
                opts.Verbose = value;

            case 'stoponfailure'
                opts.StopOnFailure = value;

            case 'erroronfailure'
                opts.ErrorOnFailure = value;

            otherwise
                error('Unknown option: %s.', name);
        end
    end
end