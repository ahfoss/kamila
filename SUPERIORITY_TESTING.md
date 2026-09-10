# Algorithm Performance Superiority Testing in `kamila`

This document details the **Performance Superiority Testing Framework** implemented in `kamila`.

---

## 1. Motivation: Why Superiority Testing?

Standard continuous integration performance testing often employs a *non-regression* design: it only checks that proposed changes do not make the code slower.

However, algorithmic optimizations often introduce additional code complexity, specialized data structures, or branch logic. Accepting pull requests that yield only negligible, sub-threshold differences (e.g., 0.5% or 1–2 milliseconds of noise) clutters the codebase and increases maintenance burden without providing practical user benefit.

To maintain a clean and performant codebase, `kamila` enforces a **Superiority Framework**:
- The **null hypothesis** assumes the candidate change is **inferior or equivalent**.
- A pull request claiming a performance optimization must **reject this null hypothesis** by demonstrating statistically significant superiority by at least a specified margin $\delta$.
- If a candidate cannot prove superiority, the PR fails the check and should not be merged.

---

## 2. Statistical Methodology

### Non-Parametric Unpaired Testing
Benchmarked execution times on modern operating systems are inherently right-skewed due to R garbage collection pauses, CPU frequency scaling, and OS process scheduling. Parametric assumptions (such as the normality required by Student's or Welch's $t$-test) are frequently violated.

We use the **two-sample Mann-Whitney U test** (Wilcoxon rank-sum test), which makes no distributional assumptions about runtime normality and evaluates stochastic dominance between candidate and baseline runtimes.

### Superiority Hypothesis with Margin $\delta$
Let $\text{Location}_{\text{base}}$ and $\text{Location}_{\text{cand}}$ denote the median runtime parameters of the baseline and candidate implementations, respectively.

For a defined superiority margin $\delta \in [0, 1)$ (default $\delta = 0.01$, or $1\%$ faster):

$$\begin{aligned}
H_0 &: \text{Location}_{\text{cand}} \ge (1 - \delta) \cdot \text{Location}_{\text{base}} \quad \text{(Inferior or Equivalent)} \\
H_1 &: \text{Location}_{\text{cand}} < (1 - \delta) \cdot \text{Location}_{\text{base}} \quad \text{(Superior by at least } \delta\text{)}
\end{aligned}$$

- **Decision Rule**: Reject $H_0$ if $p < \alpha$ (default $\alpha = 0.05$).
- When $H_0$ is rejected, the result is marked **`SUPERIOR (PASS)`**.
- If $H_0$ cannot be rejected, the result is marked **`NOT SUPERIOR (FAIL)`**.
- Per project policy, raw unadjusted $p$-values are evaluated across conditions without multiple testing corrections.

---

## 3. Dataset Conditions

Certain optimizations (such as C++ memory layout, vectorization, or algorithmic complexity improvements) only manifest at larger scales, while others may add initialization overhead that hurts small datasets. The benchmark evaluates three distinct scales:

| Tier | Sample Size ($N$) | Features ($P_{\text{con}}, P_{\text{cat}}$) | Replications | Purpose |
| :--- | :--- | :--- | :--- | :--- |
| **Small** | $500$ | $4\text{ continuous}, 4\text{ categorical}$ | **30 runs** | Verifies low overhead and consistency on small datasets. |
| **Medium** | $50,000$ | $4\text{ continuous}, 4\text{ categorical}$ | **30 runs** | Evaluates typical multi-thousand observation clustering workloads. |
| **Large** | $500,000$ | $4\text{ continuous}, 4\text{ categorical}$ | **15 runs** | Stress-tests scaling, memory bandwidth, and C++ inner loops. |

Total execution time across all 3 tiers is $\approx 4$ minutes, providing rapid feedback in CI and local development.

---

## 4. Selective CI Triggering for Pull Requests

Performance benchmarking is computationally intensive and not required for every PR (e.g., bug fixes, documentation, or statistical correctness updates).

### Triggering via PR Label
To run superiority benchmarking on a pull request:
1. Add the **`benchmark`** label to your PR in GitHub.
2. The [`.github/workflows/performance-superiority.yaml`](.github/workflows/performance-superiority.yaml) workflow will automatically launch.
3. The workflow builds `master` (baseline) and your PR branch (candidate) into isolated R library environments and executes the benchmark suite.
4. The workflow publishes a detailed Markdown report table to the GitHub Actions Job Summary.
5. If the candidate fails to demonstrate superiority, the CI check will fail, preventing the PR from being merged until the performance goal is met.

### Manual Triggering (`workflow_dispatch`)
You can also manually run the workflow against any branch or PR from the GitHub Actions tab:
1. Navigate to **Actions** $\rightarrow$ **Performance Superiority Benchmark**.
2. Click **Run workflow**, select the target branch, and optionally configure the superiority margin $\delta$ (e.g., `0.10` for a 10% threshold).

---

## 5. Running Benchmarks Locally

Before opening a pull request, you can validate your proposed optimizations on your local machine using `inst/benchmarks/run_superiority_benchmark.R`.

### A. Saving a Baseline Before Making Changes
From the clean `master` branch:
```bash
Rscript inst/benchmarks/run_superiority_benchmark.R --save-baseline baseline_master.rds
```

### B. Testing Your Optimizations Against the Baseline
Switch to your feature branch, recompile the package, and test against your saved baseline:
```bash
# Recompile package with your changes
Rscript -e "Rcpp::compileAttributes(); pkgload::load_all()"

# Run the superiority benchmark against the recorded baseline
Rscript inst/benchmarks/run_superiority_benchmark.R --baseline-file baseline_master.rds --delta 0.01
```

### C. Useful Command-Line Options

| Option | Default | Description |
| :--- | :--- | :--- |
| `--delta <num>` | `0.01` | Required superiority margin $\delta$ ($0.01 = 1\%$). |
| `--alpha <num>` | `0.05` | Significance level for hypothesis testing. |
| `--require-all <bool>` | `TRUE` | Require all executed tiers to achieve superiority to pass. |
| `--tiers <list>` | `small,medium,large` | Comma-separated list of tiers to run (e.g., `--tiers large`). |
| `--quick` | `FALSE` | Runs a small smoke test (5 small, 5 medium, 3 large runs) for quick verification. |
| `--output-md <path>` | `""` | Output path to save the generated Markdown summary table. |
