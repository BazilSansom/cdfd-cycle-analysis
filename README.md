# CDFD Cycle Analysis

MATLAB tools for Circular--Directional Flow Decomposition (CDFD), Balanced Flow Forwarding (BFF), probabilistic cycle decomposition, throughput circularity, and signed feedback analysis.

This repository accompanies ongoing research on using CDFD to separate circular and directional structure in weighted directed networks, and then to decompose the circular component into probabilistically weighted simple cycles.

## Core workflow

```matlab
% Original weighted directed network
W = [...];

% BFF CDFD decomposition
[C, D, info] = cdfd_bff(W);

% Exact probabilistic cycle decomposition of the circular component
[Tcycles, status] = cdfdExactAllCycles(C);

% Network-level summaries
[S, diagnostics] = cdfdCycleSummaries(W, C, Tcycles);
```

The high-level wrapper runs the same workflow directly:

```matlab
R = cdfdCycleAnalysis(W);
```

## Main functions

- `cdfd_bff.m`  
  Computes the BFF CDFD decomposition `W = C + D`, where `C` is balanced/circular and `D` is acyclic/directional.

- `cdfdExactAllCycles.m`  
  Enumerates all simple directed cycles in `C`, including self-loops, and computes exact Qian--Qian/Kalpazidou probabilistic cycle coefficients.

- `cdfdExactCycleWeight.m`  
  Computes the exact probabilistic coefficient for one nominated cycle.

- `cdfdCycleSummaries.m`  
  Computes summary quantities such as circular volume, circular throughput, mean cycle length, throughput circularity, edge-volume circularity, and cycle concentration.

- `cdfdCycleAnalysis.m`  
  Convenience wrapper for the full exact workflow.

- `makeCdfdExampleFigures.m`  
  Generates reproducible example figures used for validation and illustration.

  - `cdfdSignedCycleSummaries.m`  
  Computes positive and negative feedback-loop throughput and volume from signed edge labels.

- `cdfdCoParticipationKernel.m`  
  Computes the edge co-participation kernel \(K_{ef}\), showing which edges appear together in the same cycles.

- `cdfdTopCycles.m`  
  Returns a ranked table of cycles by throughput contribution, volume contribution, or length.

- `cdfdPlotCycleCoverage.m`  
  Plots cumulative top-\(k\) coverage curves for circular throughput and circular edge-volume.

- `tests/runBasicValidation.m`  
  Runs deterministic validation tests for the main toolbox functions.


## Validation

Run the basic validation suite from the repository root:

```matlab
addpath("tests")
results = runBasicValidation;
```

All tests should pass.

## Example figures

To regenerate the example figures:

```matlab
results = makeCdfdExampleFigures("figures");
```

This creates:

```text
figures/fig_bidirected_triangle.pdf
figures/fig_subdivision_invariance.pdf
figures/fig_signed_overlap.pdf
```

To plot top-cycle coverage:

```matlab
R = cdfdCycleAnalysis(W);
out = cdfdPlotCycleCoverage(R.Tcycles, ...
    'RankBy', 'throughput');
```

## Requirements

- MATLAB R2021a or later, for `allcycles`
- MATLAB graph functions including `digraph` and `conncomp`

## Status

Research code under active development. The exact all-cycle routines are intended for small networks and validation examples. Larger networks will require sampling or streaming Monte Carlo methods.

## Citation

Citation details will be added once the associated paper is available.

## Licence

Licence to be added.