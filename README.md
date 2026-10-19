# CRISPR Trace Analysis

**CRISPR Trace Analysis** is an R Shiny application for estimating CRISPR/Cas9 editing efficiency and indel distributions from Sanger sequencing chromatograms (`.ab1` files).

The application compares wild-type (WT) and CRISPR-edited sequencing traces, compensates for local differences in capillary electrophoresis, and estimates the relative contributions of different insertion and deletion sizes using non-negative least squares (NNLS) deconvolution.

CRISPR Trace Analysis estimates CRISPR editing efficiency by aligning WT and edited Sanger chromatograms and fitting WT-derived indel models to the four-channel fluorescence signals downstream of the predicted cleavage site using non-negative least squares (NNLS).

## Motivation

CRISPR Trace Analysis was developed to provide a freely available, open-source alternative for analysing Sanger sequencing traces of CRISPR-edited samples. Existing tools for estimating genome-editing efficiency often do not provide access to their source code, limiting transparency and independent verification of their underlying methods.

This project aims to make the analysis process more transparent, reproducible, and accessible by providing an openly available implementation that researchers can inspect, evaluate, adapt, and extend to suit their needs.

## Features

- **AB1 chromatogram analysis:** Load WT and CRISPR-edited Sanger sequencing files.
- **gRNA recognition:** Identify the guide RNA in the WT sequence, including reverse-complement matching, and estimate the SpCas9 cleavage site.
- **Optional reference sequence:** Supply a reference DNA sequence to validate the gRNA and recover the cleavage-site position when the WT read is truncated.
- **Flexible alignment:** Align WT and CRISPR traces using sequential peak pairing, smooth-spline warping, and pre-cut signal correlation optimization.
- **Amplitude normalization:** Normalize CRISPR trace intensities against the pre-cut region.
- **NNLS-based indel estimation:** Deconvolve the edited trace into a configurable range of insertion/deletion size models.
- **Editing efficiency estimation:** Calculate the estimated non-WT fraction from the fitted indel distribution.
- **Interactive trace visualization:** Compare WT and CRISPR traces, inspect individual indel models, and reconstruct mixtures of selected indel components.
- **Four-channel analysis:** Inspect A, C, G, and T fluorescence signals alongside the corresponding modeled traces.
- **Sequence alignment display:** Examine artificial WT-to-indel model alignments, with deletions, unknown inserted bases, and the modeled cut position highlighted.
- **Quality and alignment diagnostics:** Report pre-cut, near-cut, and post-cut correlations, alignment parameters, peak anchors, and the model's coefficient of determination (\(R^2\)).
- **Data export:** Download estimated indel frequencies as a CSV file.

## Analysis workflow

1. **Load sequencing data** — Import WT and CRISPR-edited `.ab1` files and provide the gRNA sequence.
2. **Identify the target** — Locate the gRNA, determine its orientation, and calculate the expected cleavage position. An optional reference sequence supports validation and recovery when the WT trace does not fully cover the target.
3. **Extract the analysis window** — Select a configurable region upstream and downstream of the cleavage site.
4. **Align the chromatograms** — Match sequential peak positions and fit a smooth spline to account for local differences in electrophoretic mobility. Optimize the residual offset using the pre-cut signal.
5. **Normalize signal amplitudes** — Match the average signal intensity in the pre-cut region.
6. **Model indels** — Generate WT-derived trace models with configurable downstream shifts corresponding to different indel sizes.
7. **Estimate indel frequencies** — Fit the modeled signals to the CRISPR trace using NNLS and normalize the resulting coefficients to obtain relative fractions.
8. **Inspect and export results** — Review the indel spectrum, compare four-channel traces, inspect alignment diagnostics, and export the estimated frequencies.

## Methodology

### Peak-to-peak spline alignment

The alignment procedure uses a common starting position identified through sequence alignment. Subsequent WT and CRISPR peaks are paired sequentially, assuming that corresponding peak indices represent successive sequencing products.

A smooth spline maps WT trace coordinates onto CRISPR trace coordinates, compensating for local differences in peak spacing. A residual integer offset is then optimized against the pre-cut region, where the sequence is expected to remain unchanged.

This approach aims to preserve the correspondence between sequencing products while accounting for local mobility differences. Correlations before, near, and after the cleavage site provide diagnostic measures of alignment quality.

### NNLS-based indel deconvolution

The application constructs a set of candidate chromatogram models by shifting the WT fluorescence signal downstream of the cleavage site. The shifts correspond to candidate indel sizes, including a zero-shift WT model.

Non-negative least squares is used to fit a weighted combination of these models to the aligned CRISPR trace:

$$
\mathbf{y} \approx \mathbf{X}\boldsymbol{\beta},
\qquad \boldsymbol{\beta}\geq 0
$$

where:

- \(\mathbf{y}\) is the aligned CRISPR trace, including all four fluorescence channels.
- \(\mathbf{X}\) contains the candidate WT-derived indel models.
- \(\boldsymbol{\beta}\) contains their non-negative fitted coefficients.

The coefficients are normalized to estimate the relative contribution of each modeled indel class. The editing efficiency is estimated as:

\[
\text{Editing efficiency} = 100\% - \text{WT fraction}
\]

The fit quality is summarized using \(R^2\).

**Important:** The indel models are simplified, WT-derived signal-shift models. Inserted nucleotide identities are not reconstructed, and the estimated fractions are model-dependent rather than direct measurements of individual alleles. The results should therefore be interpreted as estimates of the indel-size distribution, not as definitive sequence-level identification or quantification of every edited allele.

## Technology stack

- **R** — Numerical analysis and signal processing
- **Shiny** — Interactive web interface
- **sangerseqR** — Sanger sequencing chromatogram processing
- **Biostrings** — DNA sequence handling, alignment, and reverse-complement matching
- **nnls** — Non-negative least-squares deconvolution

## Intended use

CRISPR Trace Analysis is intended as a research tool for exploratory analysis of CRISPR/Cas9 editing outcomes from Sanger sequencing data. It provides visual and numerical diagnostics to help users assess chromatogram alignment and interpret estimated indel distributions.

Results depend on sequencing quality, target coverage, trace alignment, and the assumptions underlying the WT-derived indel models. In particular, poor-quality traces, incomplete coverage of the cleavage site, and complex editing outcomes may reduce the reliability of the estimates.

## Status

This project is under active development. The analysis methods, visualization, and diagnostic features may evolve as the implementation is refined and validated.

## License

See the repository's `LICENSE` file for licensing terms.
