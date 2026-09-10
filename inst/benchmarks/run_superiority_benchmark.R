#!/usr/bin/env Rscript

# ==============================================================================
# run_superiority_benchmark.R
#
# Performance Superiority Testing Suite for kamila
# Evaluates whether a candidate version is statistically significantly superior
# (faster) than a baseline version using an unpaired non-parametric Mann-Whitney
# U test with a superiority margin delta:
#
#   H0: Location(Candidate) >= (1 - delta) * Location(Baseline)
#   H1: Location(Candidate) <  (1 - delta) * Location(Baseline)
#
# Dataset Conditions:
#   - Small:  N =     500, 30 runs
#   - Medium: N =  50,000, 30 runs
#   - Large:  N = 500,000, 15 runs
# ==============================================================================

# ------------------------------------------------------------------------------
# Tier Definitions
# ------------------------------------------------------------------------------
TIER_CONFIGS <- list(
  small = list(
    name = "Small",
    sampSize = 500,
    nConVar = 4,
    nCatVar = 4,
    nCatLevels = 4,
    numClust = 2,
    numInit = 5,
    maxIter = 15,
    runs = 30,
    quick_runs = 5
  ),
  medium = list(
    name = "Medium",
    sampSize = 50000,
    nConVar = 4,
    nCatVar = 4,
    nCatLevels = 4,
    numClust = 2,
    numInit = 1,
    maxIter = 15,
    runs = 30,
    quick_runs = 5
  ),
  large = list(
    name = "Large",
    sampSize = 500000,
    nConVar = 4,
    nCatVar = 4,
    nCatLevels = 4,
    numClust = 2,
    numInit = 1,
    maxIter = 10,
    runs = 15,
    quick_runs = 3
  )
)

# ------------------------------------------------------------------------------
# Worker Mode Execution
# Runs within an isolated Rscript sub-process to eliminate namespace/DLL locks
# ------------------------------------------------------------------------------
run_worker <- function(lib_path, tier_name, runs, seed) {
  if (!is.na(lib_path) && nzchar(lib_path)) {
    .libPaths(c(normalizePath(lib_path, mustWork = TRUE), .libPaths()))
  }

  if (!requireNamespace("kamila", quietly = TRUE)) {
    stop("Package 'kamila' could not be loaded from library: ", lib_path)
  }

  tier <- TIER_CONFIGS[[tier_name]]
  if (is.null(tier)) {
    stop("Unknown tier: ", tier_name)
  }

  set.seed(seed)
  runtimes <- numeric(runs)

  for (i in seq_len(runs)) {
    # Generate mixed data
    dat <- kamila::genMixedData(
      sampSize = tier$sampSize,
      nConVar = tier$nConVar,
      nCatVar = tier$nCatVar,
      nCatLevels = tier$nCatLevels,
      nConWithErr = 1,
      nCatWithErr = 1,
      popProportions = c(0.5, 0.5),
      conErrLev = 0.2,
      catErrLev = 0.2
    )

    conDf <- as.data.frame(dat$conVars)
    catDf <- as.data.frame(lapply(as.data.frame(dat$catVars), factor))

    # Clean garbage collection prior to timing
    gc(verbose = FALSE)

    t0 <- proc.time()
    res <- kamila::kamila(
      conVar = conDf,
      catFactor = catDf,
      numClust = tier$numClust,
      numInit = tier$numInit,
      maxIter = tier$maxIter,
      verbose = FALSE,
      calcNumClust = "none"
    )
    t_elapsed <- (proc.time() - t0)[3]
    runtimes[i] <- t_elapsed
  }

  # Output comma-separated runtimes
  cat(paste(sprintf("%.6f", runtimes), collapse = ","))
}

# ------------------------------------------------------------------------------
# Coordinator Execution Helpers
# ------------------------------------------------------------------------------
invoke_worker <- function(rscript, script_path, lib_path, tier_name, runs, seed) {
  args <- c(
    "--vanilla",
    script_path,
    "--worker",
    "--tier", tier_name,
    "--runs", as.character(runs),
    "--seed", as.character(seed)
  )
  if (!is.na(lib_path) && nzchar(lib_path)) {
    args <- c(args, "--lib", normalizePath(lib_path, mustWork = TRUE))
  }

  out <- system2(rscript, args = args, stdout = TRUE, stderr = TRUE)
  out_lines <- out[nzchar(trimws(out))]
  if (length(out_lines) == 0) {
    stop("Worker produced no output. Output was:\n", paste(out, collapse = "\n"))
  }

  times_str <- out_lines[length(out_lines)]
  times <- as.numeric(strsplit(times_str, ",")[[1]])
  if (any(is.na(times)) || length(times) != runs) {
    stop("Failed to parse runtimes from worker output:\n", paste(out, collapse = "\n"))
  }
  times
}

# ------------------------------------------------------------------------------
# Statistical Analysis: Mann-Whitney U Superiority Test with Margin delta
# ------------------------------------------------------------------------------
analyze_tier <- function(cand_times, base_times, delta = 0.01, alpha = 0.05) {
  # Null hypothesis: Location(cand) >= (1 - delta) * Location(base)
  # Alternative:     Location(cand) <  (1 - delta) * Location(base) (superior)
  scaled_base <- base_times * (1 - delta)

  wt <- wilcox.test(cand_times, scaled_base, alternative = "less", exact = FALSE)

  med_base <- median(base_times)
  iqr_base <- IQR(base_times)
  mean_base <- mean(base_times)
  sd_base <- sd(base_times)

  med_cand <- median(cand_times)
  iqr_cand <- IQR(cand_times)
  mean_cand <- mean(cand_times)
  sd_cand <- sd(cand_times)

  pct_speedup <- (med_base - med_cand) / med_base * 100
  mean_pct_speedup <- (mean_base - mean_cand) / mean_base * 100

  p_val <- wt$p.value
  is_superior <- (p_val < alpha)

  list(
    n_base = length(base_times),
    n_cand = length(cand_times),
    med_base = med_base,
    iqr_base = iqr_base,
    mean_base = mean_base,
    sd_base = sd_base,
    med_cand = med_cand,
    iqr_cand = iqr_cand,
    mean_cand = mean_cand,
    sd_cand = sd_cand,
    pct_speedup = pct_speedup,
    mean_pct_speedup = mean_pct_speedup,
    statistic = as.numeric(wt$statistic),
    p_value = p_val,
    is_superior = is_superior
  )
}

format_time <- function(val_sec) {
  if (val_sec < 1) {
    sprintf("%.1f ms", val_sec * 1000)
  } else {
    sprintf("%.2f s", val_sec)
  }
}

# ------------------------------------------------------------------------------
# Main Coordinator Routine
# ------------------------------------------------------------------------------
main <- function() {
  args <- commandArgs(trailingOnly = TRUE)

  # Parse arguments
  parsed <- list(
    worker = FALSE,
    lib = "",
    tier = "",
    runs = 0,
    seed = 42,
    base_lib = "",
    cand_lib = "",
    baseline_file = "",
    save_baseline = "",
    delta = 0.01,
    alpha = 0.05,
    require_all = TRUE,
    quick = FALSE,
    tiers = "small,medium,large",
    output_md = "",
    verbose = FALSE
  )

  i <- 1
  while (i <= length(args)) {
    arg <- args[i]
    if (arg == "--worker") {
      parsed$worker <- TRUE
    } else if (arg == "--lib" && i < length(args)) {
      parsed$lib <- args[i + 1]
      i <- i + 1
    } else if (arg == "--tier" && i < length(args)) {
      parsed$tier <- args[i + 1]
      i <- i + 1
    } else if (arg == "--runs" && i < length(args)) {
      parsed$runs <- as.integer(args[i + 1])
      i <- i + 1
    } else if (arg == "--seed" && i < length(args)) {
      parsed$seed <- as.integer(args[i + 1])
      i <- i + 1
    } else if (arg == "--base-lib" && i < length(args)) {
      parsed$base_lib <- args[i + 1]
      i <- i + 1
    } else if (arg == "--cand-lib" && i < length(args)) {
      parsed$cand_lib <- args[i + 1]
      i <- i + 1
    } else if (arg == "--baseline-file" && i < length(args)) {
      parsed$baseline_file <- args[i + 1]
      i <- i + 1
    } else if (arg == "--save-baseline" && i < length(args)) {
      parsed$save_baseline <- args[i + 1]
      i <- i + 1
    } else if (arg == "--delta" && i < length(args)) {
      parsed$delta <- as.numeric(args[i + 1])
      i <- i + 1
    } else if (arg == "--alpha" && i < length(args)) {
      parsed$alpha <- as.numeric(args[i + 1])
      i <- i + 1
    } else if (arg == "--require-all" && i < length(args)) {
      parsed$require_all <- as.logical(args[i + 1])
      i <- i + 1
    } else if (arg == "--quick") {
      parsed$quick <- TRUE
    } else if (arg == "--tiers" && i < length(args)) {
      parsed$tiers <- args[i + 1]
      i <- i + 1
    } else if (arg == "--output-md" && i < length(args)) {
      parsed$output_md <- args[i + 1]
      i <- i + 1
    } else if (arg == "--verbose") {
      parsed$verbose <- TRUE
    }
    i <- i + 1
  }

  # If running in worker mode, execute and exit immediately
  if (parsed$worker) {
    run_worker(
      lib_path = parsed$lib,
      tier_name = parsed$tier,
      runs = parsed$runs,
      seed = parsed$seed
    )
    return(invisible(0))
  }

  # --- Coordinator Mode ---
  cat("======================================================================\n")
  cat("KAMILA Algorithm Performance Superiority Testing Suite\n")
  cat("======================================================================\n")
  cat(sprintf("Superiority Margin (delta) : %.1f%%\n", parsed$delta * 100))
  cat(sprintf("Significance Level (alpha) : %.3f\n", parsed$alpha))
  cat(sprintf("Hypothesis                 : H0: Location(Cand) >= %.3f * Location(Base)\n", 1 - parsed$delta))
  cat(sprintf(
    "                             H1: Location(Cand) <  %.3f * Location(Base) [Superior]\n",
    1 - parsed$delta
  ))
  cat(sprintf("Require All Tiers Superior : %s\n", parsed$require_all))
  if (parsed$quick) cat("Mode                       : Quick Smoke Test\n")
  cat("----------------------------------------------------------------------\n\n")

  active_tiers <- strsplit(parsed$tiers, "[, ]+")[[1]]
  active_tiers <- active_tiers[nzchar(active_tiers)]

  # Identify Rscript and this script's path
  rscript <- file.path(R.home("bin"), "Rscript")
  script_args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", script_args, value = TRUE)
  if (length(file_arg) > 0) {
    script_path <- sub("^--file=", "", file_arg[1])
  } else {
    script_path <- "inst/benchmarks/run_superiority_benchmark.R"
  }

  # Baseline timings storage
  baseline_results <- list()
  if (nzchar(parsed$baseline_file) && file.exists(parsed$baseline_file)) {
    cat(sprintf("Loading saved baseline timings from: %s\n", parsed$baseline_file))
    baseline_results <- readRDS(parsed$baseline_file)
  }

  results <- list()
  all_passed <- TRUE
  any_passed <- FALSE

  for (t_name in active_tiers) {
    tier <- TIER_CONFIGS[[t_name]]
    if (is.null(tier)) {
      warning("Unknown tier: ", t_name, ", skipping.")
      next
    }

    total_runs <- if (parsed$quick) tier$quick_runs else tier$runs
    cat(sprintf("Executing Tier: %s (N = %s, Target Runs = %d)...\n",
                tier$name, format(tier$sampSize, big.mark = ","), total_runs))

    # 1. Obtain Baseline Timings
    base_times <- NULL
    if (!is.null(baseline_results[[t_name]])) {
      base_times <- baseline_results[[t_name]]
      cat(sprintf("  -> Using %d pre-recorded baseline runs.\n", length(base_times)))
    } else {
      cat("  -> Running baseline variant...")
      base_times <- invoke_worker(
        rscript = rscript,
        script_path = script_path,
        lib_path = parsed$base_lib,
        tier_name = t_name,
        runs = total_runs,
        seed = parsed$seed
      )
      cat(sprintf(" Done (Median = %s).\n", format_time(median(base_times))))
      baseline_results[[t_name]] <- base_times
    }

    # If only saving baseline, skip candidate
    if (nzchar(parsed$save_baseline) && !nzchar(parsed$cand_lib)) {
      next
    }

    # 2. Obtain Candidate Timings
    cat("  -> Running candidate variant...")
    cand_times <- invoke_worker(
      rscript = rscript,
      script_path = script_path,
      lib_path = parsed$cand_lib,
      tier_name = t_name,
      runs = total_runs,
      seed = parsed$seed + 1000
    )
    cat(sprintf(" Done (Median = %s).\n", format_time(median(cand_times))))

    # 3. Statistical Analysis
    analysis <- analyze_tier(cand_times, base_times, delta = parsed$delta, alpha = parsed$alpha)
    analysis$tier_name <- tier$name
    analysis$sampSize <- tier$sampSize
    results[[t_name]] <- analysis

    if (analysis$is_superior) {
      any_passed <- TRUE
      cat(sprintf("  -> Result: SUPERIOR (Speedup: %+.1f%%, p-val = %.4e)\n\n",
                  analysis$pct_speedup, analysis$p_value))
    } else {
      all_passed <- FALSE
      cat(sprintf("  -> Result: NOT SUPERIOR (Speedup: %+.1f%%, p-val = %.4f)\n\n",
                  analysis$pct_speedup, analysis$p_value))
    }
  }

  # Save baseline if requested
  if (nzchar(parsed$save_baseline)) {
    cat(sprintf("Saving baseline measurements to: %s\n", parsed$save_baseline))
    saveRDS(baseline_results, parsed$save_baseline)
    if (!nzchar(parsed$cand_lib)) {
      cat("Baseline collection complete.\n")
      return(invisible(0))
    }
  }

  # ----------------------------------------------------------------------------
  # Generate Markdown Summary
  # ----------------------------------------------------------------------------
  final_verdict <- if (parsed$require_all) all_passed else any_passed
  verdict_str <- if (final_verdict) "PASSED - SUPERIORITY DEMONSTRATED" else "FAILED - INFERIOR OR EQUIVALENT"

  verdict_badge <- if (final_verdict) "**PASS** :white_check_mark:" else "**FAIL** :x:"
  test_desc <- paste0(
    "- **Statistical Test:** Unpaired Mann-Whitney U test ",
    "($H_0: \\text{Location}_{\\text{cand}} \\ge (1 - \\delta) \\cdot \\text{Location}_{\\text{base}}$)"
  )
  table_header <- paste(
    "| Condition | Sample Size ($N$) | Runs | Baseline Median (IQR) |",
    "Candidate Median (IQR) | Median Speedup | $p$-value | Verdict |"
  )
  table_sep <- paste(
    "| :--- | :--- | :--- | :--- |",
    ":--- | :--- | :--- | :--- |"
  )

  md <- c(
    "## Algorithm Performance Superiority Test Results",
    "",
    sprintf("- **Overall Verdict:** %s", verdict_badge),
    sprintf("- **Superiority Margin ($\\delta$):** %.1f%% faster than baseline required", parsed$delta * 100),
    test_desc,
    sprintf("- **Significance Threshold ($\\alpha$):** %.3f (no multiplicity adjustment)", parsed$alpha),
    "",
    table_header,
    table_sep
  )

  for (t_name in names(results)) {
    res <- results[[t_name]]
    p_str <- if (res$p_value < 1e-4) sprintf("%.2e", res$p_value) else sprintf("%.4f", res$p_value)
    status_icon <- if (res$is_superior) "**SUPERIOR** :white_check_mark:" else "**NOT SUPERIOR** :x:"
    speedup_str <- sprintf("%+.1f%%", res$pct_speedup)
    if (res$is_superior) speedup_str <- paste0("**", speedup_str, "**")

    row <- paste(
      sprintf(
        "| **%s** | %s | %d / %d |",
        res$tier_name,
        format(res$sampSize, big.mark = ","),
        res$n_base,
        res$n_cand
      ),
      sprintf(
        "%s (%s) | %s (%s) |",
        format_time(res$med_base),
        format_time(res$iqr_base),
        format_time(res$med_cand),
        format_time(res$iqr_cand)
      ),
      sprintf("%s | `%s` | %s |", speedup_str, p_str, status_icon)
    )
    md <- c(md, row)
  }
  md <- c(md, "")

  md_text <- paste(md, collapse = "\n")

  # Print to stdout
  cat("----------------------------------------------------------------------\n")
  cat(md_text, "\n")
  cat("======================================================================\n")
  cat(sprintf("Final Verdict: %s\n", verdict_str))
  cat("======================================================================\n")

  # Write to file if requested
  if (nzchar(parsed$output_md)) {
    writeLines(md_text, parsed$output_md)
    cat(sprintf("Wrote report to: %s\n", parsed$output_md))
  }

  # Write to GitHub Step Summary if running in GitHub Actions
  step_summary <- Sys.getenv("GITHUB_STEP_SUMMARY", unset = "")
  if (nzchar(step_summary)) {
    write(md_text, file = step_summary, append = TRUE)
  }

  # Exit status: 0 if passed, 1 if failed
  if (!final_verdict) {
    quit(status = 1, save = "no")
  }
}

if (!interactive()) {
  main()
}
