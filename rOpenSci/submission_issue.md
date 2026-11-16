Submitting Author Name: Jean Pauphilet
Submitting Author Github Handle: <!--author1-->@jeanpauphilet<!--end-author1-->
Other Package Authors Github handles: (comma separated, delete if none) <!--author-others-->[@ryan_github_handle]<!--end-author-others-->
Repository:  <!--repourl-->https://github.com/jeanpauphilet/msPCA<!--end-repourl-->
Version submitted: 0.5.1
Submission type: <!--submission-type-->Stats<!--end-submission-type-->
Badge grade: <!--statsgrade-->silver<!--end-statsgrade-->
Editor: <!--editor--> TBD <!--end-editor-->
Reviewers: <!--reviewers-list--> TBD <!--end-reviewers-list-->
<!--due-dates-list--><!--end-due-dates-list-->
Archive: TBD
Version accepted: TBD
Language: <!--language-->en<!--end-language-->

---



-   Paste the full DESCRIPTION file inside a code block below:

```
[paste the output of: cat DESCRIPTION]
```

## Scope 

- Please indicate which of our [statistical package categories](https://stats-devguide.ropensci.org/overview.html#overview-categories) this package falls under:

     **Statistical Packages**

	- [ ] Bayesian and Monte Carlo Routines
	- [x] Dimensionality Reduction, Clustering, and Unsupervised Learning
	- [ ] Machine Learning
	- [ ] Regression and Supervised Learning
	- [ ] Exploratory Data Analysis (EDA) and Summary Statistics
	- [ ] Spatial Analyses
	- [ ] Time Series Analyses
	- [ ] Probability Distributions

## Pre-submission Inquiry

- [ ] A pre-submission inquiry has been approved in [issue#<issue_num>](https://github.com/ropensci/software-review/issues/<issue_num>)
  [None made; delete this section if so.]


## General Information


-   Who is the target audience and what are scientific applications of this package?

    Researchers and practitioners who need interpretable dimensionality reduction or unsupervised feature selection on moderate to high-dimensional data. Examples include genomics, finance (the replication study uses S&P 500 returns) and the social sciences. msPCA computes several sparse principal components jointly. Each component has a prescribed number of non-zero loadings, and the components are (near-)orthogonal.

-   Paste your responses to our [*General Standard* **G1.1**](https://stats-devguide.ropensci.org/standards.html#general-standards) here:

    msPCA is **the first implementation of a novel algorithm**. It implements the method of Cory-Wright and Pauphilet, "Sparse PCA with Multiple Components", *Operations Research* (2026), [doi:10.1287/opre.2023.0598](https://doi.org/10.1287/opre.2023.0598). The method uses an iterative deflation heuristic with a truncated power method at each step.

    Related software, all compared against msPCA in the [accompanying paper (arXiv:2607.05229)](https://arxiv.org/abs/2607.05229) and in the [benchmarking article](https://jeanpauphilet.github.io/msPCA/articles/benchmarking.html):
    - R: [`elasticnet`](https://cran.r-project.org/package=elasticnet) (Zou et al., 2006; elastic-net regularization, LARS-like)
    - R: [`PMA`](https://cran.r-project.org/package=PMA) (Witten et al., 2009; penalized matrix decomposition, sequential deflation)
    - R: [`sparsepca`](https://cran.r-project.org/package=sparsepca) (Erichson et al., 2020; variable projection with an l1 penalty)
    - R: [`mixOmics`](https://cran.r-project.org/package=mixOmics) (Rohart et al., 2017; regularized SVD after Shen and Huang, 2008)
    - R: [`nsprcomp`](https://cran.r-project.org/package=nsprcomp) (Sigg, 2019; EM with deflation, or a cumulative formulation)
    - R: [`amanpg`](https://cran.r-project.org/package=amanpg) (Chen et al., 2020; alternating manifold proximal gradient)
    - Python: [`scikit-learn` `SparsePCA`](https://scikit-learn.org/stable/modules/generated/sklearn.decomposition.SparsePCA.html) (Pedregosa et al., 2011; dictionary learning)

    msPCA sets an exact cardinality on each component and is the only one of these packages that supports both orthogonal loadings and uncorrelated PCs as non-redundancy constraints. It can also work directly on the data matrix, without forming the covariance matrix.

-   Versioning: version 0.5.1 is submitted. The post-review release will be 0.6.0, because the development version contains small API changes made while preparing for review: logical flags are now strictly logical, `type` and `divisor` are case-insensitive, `mspca()`/`tpm()` accept all-numeric data.frames, and results include a `control` field. See `NEWS.md`.

-   Ethics, data privacy and human subjects research: not applicable.

## Badging

-    What grade of badge are you aiming for? **Silver**

-    Which of the four aspects for silver does the package fulfil?

     - **Compliance with a good number of standards beyond the minimum.** 72 standards are documented as complied with and 30 as not applicable, per `srr_stats_pre_submit()`. [Check against the minimal vs. full applicable set.]
     - **Excellence in multiple sub-categories.** [Optional: documentation (README, vignette, replication study, pkgdown site) and testing (about 98% coverage, correctness checks against exhaustive search and eigen-decomposition, parameter-recovery and noise tests).]
     - **Generality of usage.** [Optional: both covariance (`type = "Sigma"`) and data (`type = "X"`) input, any correlation estimator, optional feasibility constraints.]
     - **Package structure and design.** [Leave to reviewers.]

     Deliberately not met (each has an `@srrstatsNA` justification in `R/srr-stats-standards.R`): UL6.0 (no `plot()` method), UL3.3 and UL7.4 (no `predict()` method; scores are `newX %*% x_best`), G5.10–G5.12 (no extended tests), and G5.4b/G5.4c (the method is novel, so there is no prior implementation to compare against).


## Technical checks

- [x] I have read the [rOpenSci packaging guide](https://devguide.ropensci.org/building.html).
- [x] I have read the [author guide](https://devguide.ropensci.org/softwarereview_author.html) and I expect to maintain this package for at least 2 years or have another maintainer identified.
- [x] I/we have read the [*Statistical Software Peer Review* Guide for Authors](https://stats-devguide.ropensci.org/pkgdev.html).
- [x] I/we have run [`autotest`](https://github.com/ropensci-review-tools/autotest) checks on the package, and ensured no tests fail. [Re-run on the final commit.]
- [x] The [`srr_stats_pre_submit()` function](https://ropensci-review-tools.github.io/srr/reference/srr_stats_pre_submit.html) confirms this package may be submitted.
- [x] The [`pkgcheck()` function](https://docs.ropensci.org/pkgcheck/reference/pkgcheck.html) confirms this package may be submitted. [Attach or paste the output.]

This package:

- [x] does not violate the Terms of Service of any service it interacts with.
- [x] has a CRAN and OSI accepted license (MIT).
- [x] contains a [README with instructions for installing the development version](https://devguide.ropensci.org/building.html#readme).

## Use of Generative AI

- [x] Generative AI tools were used to produce some of the material in this submission.

GitHub Copilot (an AI coding assistant) was used to help with:
- input-validation hardening after `autotest`;
- tests in `tests/testthat/test-standards.R`;
- the `srr` standards annotations in `R/srr-stats-standards.R`;
- documentation and NEWS updates.

All AI-assisted changes were reviewed by the authors. [Confirm this is accurate and link the relevant commits.]

## Publication options

- [x] Do you intend for this package to go on CRAN? (It is already on CRAN: msPCA 0.5.1.)
- [ ] Do you intend for this package to go on Bioconductor?


## Code of conduct

- [x] I agree to abide by [rOpenSci's Code of Conduct](https://ropensci.org/code-of-conduct/) during the review process and in maintaining my package should it be accepted.
