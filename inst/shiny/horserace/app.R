# ==============================================================================
# Mixed-Type Clustering Horse-Race - Interactive Shiny Web Application
# Compatible with local Shiny and Shinylive / WebR (Client-Side Browser Execution)
# ==============================================================================

if (requireNamespace("shiny", quietly = TRUE)) {
  library(shiny)
}

# If running inside WebR / WebAssembly browser environment, automatically install packages
if (exists("webr", envir = .GlobalEnv) || Sys.getenv("WEBR") == "1") {
  tryCatch({
    webr::install(c("clustMixType", "clustrd", "VarSelLCM", "flexmix", "mvtnorm", "mclust", "MASS"))
  }, error = function(e) NULL)
}

# Optional dependencies loaded conditionally for WebR compatibility
has_pkg <- function(pkg) requireNamespace(pkg, quietly = TRUE)

# ------------------------------------------------------------------------------
# Helpers: Safe Scaling & Evaluation Metrics
# ------------------------------------------------------------------------------
safe_scale <- function(m) {
  m_mat <- as.matrix(m)
  sds <- apply(m_mat, 2, sd, na.rm = TRUE)
  sds_valid <- ifelse(is.na(sds) | sds == 0, 1, sds)
  m_scaled <- scale(m_mat, center = TRUE, scale = sds_valid)
  m_scaled[is.na(m_scaled)] <- 0
  m_scaled
}

calc_ari <- function(true_labels, pred_labels) {
  true_v <- as.vector(as.integer(true_labels))
  pred_v <- as.vector(as.integer(pred_labels))

  if (has_pkg("mclust")) {
    return(mclust::adjustedRandIndex(true_v, pred_v))
  }
  tab <- table(true_v, pred_v)
  n <- length(true_v)
  if (n < 2) return(1.0)

  comb2 <- function(x) x * (x - 1) / 2
  sum_comb_tab <- sum(comb2(tab))
  sum_comb_rows <- sum(comb2(rowSums(tab)))
  sum_comb_cols <- sum(comb2(colSums(tab)))

  expected <- (sum_comb_rows * sum_comb_cols) / comb2(n)
  max_val <- 0.5 * (sum_comb_rows + sum_comb_cols)

  if (max_val == expected) return(1.0)
  (sum_comb_tab - expected) / (max_val - expected)
}

calc_misclass_error <- function(true_labels, pred_labels) {
  true_v <- as.vector(as.integer(true_labels))
  pred_v <- as.vector(as.integer(pred_labels))

  tab <- table(true_v, pred_v)
  k_true <- nrow(tab)
  k_pred <- ncol(tab)

  if (k_pred > 6 || k_true > 6) {
    # Greedy matching heuristic for larger K
    matched_correct <- 0
    temp_tab <- tab
    for (i in 1:min(k_true, k_pred)) {
      max_idx <- which(temp_tab == max(temp_tab), arr.ind = TRUE)[1, ]
      matched_correct <- matched_correct + temp_tab[max_idx[1], max_idx[2]]
      temp_tab[max_idx[1], ] <- -1
      temp_tab[, max_idx[2]] <- -1
    }
    return(1 - (matched_correct / length(true_v)))
  }

  # Exact permutation matching for K <= 6
  perms <- function(v) {
    if (length(v) <= 1) return(matrix(v, 1, 1))
    do.call(rbind, lapply(seq_along(v), function(i) {
      cbind(v[i], perms(v[-i]))
    }))
  }

  all_p <- perms(1:k_true)
  max_correct <- 0
  for (i in seq_len(nrow(all_p))) {
    p <- all_p[i, ]
    cols_to_use <- p[seq_len(min(k_true, k_pred))]
    cur_correct <- sum(sapply(seq_along(cols_to_use), function(idx) {
      if (idx <= ncol(tab)) tab[cols_to_use[idx], idx] else 0
    }))
    if (cur_correct > max_correct) max_correct <- cur_correct
  }
  1 - (max_correct / length(true_v))
}

# ------------------------------------------------------------------------------
# Synthetic Mixed Data Generator (1,000 to 100,000 observations)
# ------------------------------------------------------------------------------
generate_synthetic_mixed_data <- function(n = 1000, p_con = 10, p_cat = 10, k = 4,
                                          separation = 2.0, num_levels = 10, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)

  props <- rep(1 / k, k)
  cluster_assign <- sample(seq_len(k), size = n, replace = TRUE, prob = props)

  # Continuous variables: Gaussian cluster centers with controlled separation
  con_data <- matrix(0, nrow = n, ncol = p_con)
  colnames(con_data) <- paste0("Con_", seq_len(p_con))

  for (j in seq_len(p_con)) {
    centers <- seq(-separation * (k - 1) / 2, separation * (k - 1) / 2, length.out = k)
    con_data[, j] <- rnorm(n, mean = centers[cluster_assign], sd = 1.0)
  }

  # Categorical variables: Multinomial with dominant level per cluster
  cat_list <- list()
  for (j in seq_len(p_cat)) {
    col_vals <- integer(n)
    for (clust in seq_len(k)) {
      idx <- which(cluster_assign == clust)
      if (length(idx) > 0) {
        prob_vec <- rep(0.1, num_levels)
        dom_level <- ((clust - 1 + j - 1) %% num_levels) + 1
        prob_vec[dom_level] <- 0.6
        prob_vec <- prob_vec / sum(prob_vec)
        col_vals[idx] <- sample(seq_len(num_levels), size = length(idx), replace = TRUE, prob = prob_vec)
      }
    }
    cat_list[[paste0("Cat_", j)]] <- factor(paste0("L", col_vals), levels = paste0("L", seq_len(num_levels)))
  }
  cat_df <- as.data.frame(cat_list)

  list(
    conVars = as.data.frame(con_data),
    catVars = cat_df,
    trueID = cluster_assign,
    fullData = cbind(as.data.frame(con_data), cat_df)
  )
}

# ------------------------------------------------------------------------------
# Method Metadata & Execution Runners
# ------------------------------------------------------------------------------
method_meta <- list(
  kamila = list(name = "KAMILA", pkg = "kamila"),
  gower_pam = list(name = "Gower + PAM", pkg = "cluster"),
  kproto = list(name = "K-Prototypes", pkg = "clustMixType"),
  varsellcm = list(name = "VarSelLCM", pkg = "VarSelLCM"),
  flexmix_multinom = list(name = "FlexMix (Multinomial)", pkg = "flexmix"),
  flexmix_binary = list(name = "FlexMix (Binary Levels)", pkg = "flexmix")
)

detected_cores_count <- tryCatch({
  p_cores <- parallel::detectCores(logical = FALSE)
  if (is.na(p_cores) || p_cores < 1) {
    p_cores <- parallel::detectCores()
  }
  if (is.na(p_cores) || p_cores < 1) 2 else p_cores
}, error = function(e) 2)

run_single_fixed_k <- function(m, dat, k) {
  m_info <- method_meta[[m]]
  m_name <- m_info$name
  m_pkg <- m_info$pkg

  if (!has_pkg(m_pkg)) {
    return(data.frame(
      Method = m_name,
      Package = m_pkg,
      Time_ms = NA,
      ARI = NA,
      Error_Rate = NA,
      Status = "Package Not Installed",
      stringsAsFactors = FALSE
    ))
  }

  t_start <- proc.time()
  tryCatch({
    if (m == "kamila") {
      res <- kamila::kamila(
        dat$conVars, dat$catVars, numClust = k, numInit = 5, maxIter = 25, calcNumClust = "none"
      )
      memb <- as.integer(res$finalMemb)
    } else if (m == "gower_pam") {
      gdist <- cluster::daisy(dat$fullData, metric = "gower")
      pam_fit <- cluster::pam(gdist, k = k, diss = TRUE)
      memb <- as.integer(pam_fit$clustering)
    } else if (m == "kproto") {
      kp_fit <- clustMixType::kproto(dat$fullData, k = k, nstart = 3, verbose = FALSE)
      memb <- as.integer(kp_fit$cluster)
    } else if (m == "varsellcm") {
      v_fit <- VarSelLCM::VarSelCluster(
        x = dat$fullData, gvals = k, vbleSelec = FALSE, crit.varsel = "BIC", nbcores = 1
      )
      memb <- as.integer(VarSelLCM::fitted(v_fit, type = "partition"))
    } else if (m == "flexmix_multinom") {
      con_cols <- colnames(dat$conVars)
      cat_cols <- colnames(dat$catVars)
      con_form <- stats::as.formula(paste0("cbind(", paste(con_cols, collapse = ", "), ") ~ 1"))
      con_mod <- flexmix::FLXMCmvnorm(con_form, diagonal = TRUE)
      cat_mods <- lapply(cat_cols, function(col) {
        flexmix::FLXMRmultinom(stats::as.formula(paste0(col, " ~ 1")))
      })
      f_fit <- flexmix::flexmix(
        stats::as.formula(paste0(con_cols[1], " ~ 1")),
        data = dat$fullData,
        k = k,
        model = c(list(con_mod), cat_mods),
        control = list(iter.max = 25, minprior = 0.05, verbose = 0)
      )
      memb <- as.integer(flexmix::clusters(f_fit))
    } else if (m == "flexmix_binary") {
      con_mat <- safe_scale(dat$conVars)
      cat_dummy <- stats::model.matrix(~ . - 1, data = dat$catVars)
      df_comb <- as.data.frame(cbind(con_mat, cat_dummy))
      con_cols <- colnames(con_mat)
      bin_cols <- colnames(cat_dummy)
      con_form <- stats::as.formula(paste0("cbind(", paste(con_cols, collapse = ", "), ") ~ 1"))
      bin_form <- stats::as.formula(paste0("cbind(", paste(bin_cols, collapse = ", "), ") ~ 1"))
      f_fit <- flexmix::flexmix(
        stats::as.formula(paste0(con_cols[1], " ~ 1")),
        data = df_comb,
        k = k,
        model = list(
          flexmix::FLXMCmvnorm(con_form, diagonal = TRUE),
          flexmix::FLXMCmvbinary(bin_form)
        ),
        control = list(iter.max = 25, minprior = 0.05, verbose = 0)
      )
      memb <- as.integer(flexmix::clusters(f_fit))
    }

    t_elapsed <- (proc.time() - t_start)[["elapsed"]] * 1000
    ari <- calc_ari(dat$trueID, memb)
    err <- calc_misclass_error(dat$trueID, memb)

    data.frame(
      Method = m_name,
      Package = m_pkg,
      Time_ms = round(t_elapsed, 1),
      ARI = round(ari, 4),
      Error_Rate = round(err, 4),
      Status = "Success",
      stringsAsFactors = FALSE
    )
  }, error = function(e) {
    data.frame(
      Method = m_name,
      Package = m_pkg,
      Time_ms = NA,
      ARI = NA,
      Error_Rate = NA,
      Status = paste("Error:", substr(e$message, 1, 28)),
      stringsAsFactors = FALSE
    )
  })
}

run_single_selection <- function(m, dat, true_k, ps_cores = 1) {
  m_info <- method_meta[[m]]
  m_name <- m_info$name
  m_pkg <- m_info$pkg

  if (!has_pkg(m_pkg)) return(NULL)

  k_max_search <- min(10, max(5, true_k + 2))
  k_range <- 2:k_max_search

  t_start <- proc.time()
  tryCatch({
    if (m == "kamila") {
      kam_args <- list(
        conVar = dat$conVars,
        catFactor = dat$catVars,
        numClust = k_range,
        numInit = 3,
        calcNumClust = "ps",
        numPredStrCvRun = 5,
        predStrThresh = 0.6
      )
      if ("numCores" %in% names(formals(kamila::kamila))) {
        num_cores_ps <- if (!is.null(ps_cores) && ps_cores > 1) as.integer(ps_cores) else 1
        kam_args$numCores <- num_cores_ps
      }
      res_k <- do.call(kamila::kamila, kam_args)
      t_elapsed <- (proc.time() - t_start)[["elapsed"]] * 1000
      pred_k <- if (is.list(res_k$nClust)) res_k$nClust$bestNClust else res_k$nClust
      ari <- calc_ari(dat$trueID, res_k$finalMemb)
      data.frame(
        Method = m_name,
        Package = m_pkg,
        True_K = true_k,
        Predicted_K = pred_k,
        Criterion = "Prediction Strength",
        ARI = round(ari, 4),
        Time_ms = round(t_elapsed, 1),
        stringsAsFactors = FALSE
      )
    } else if (m == "gower_pam") {
      g_dist <- cluster::daisy(dat$fullData, metric = "gower")
      sils <- sapply(k_range, function(ki) {
        cluster::pam(g_dist, k = ki, diss = TRUE)$silinfo$avg.width
      })
      best_ki <- k_range[which.max(sils)]
      fit_best <- cluster::pam(g_dist, k = best_ki, diss = TRUE)
      t_elapsed <- (proc.time() - t_start)[["elapsed"]] * 1000
      ari <- calc_ari(dat$trueID, as.integer(fit_best$clustering))
      data.frame(
        Method = m_name,
        Package = m_pkg,
        True_K = true_k,
        Predicted_K = best_ki,
        Criterion = "Avg Silhouette Width",
        ARI = round(ari, 4),
        Time_ms = round(t_elapsed, 1),
        stringsAsFactors = FALSE
      )
    } else if (m == "kproto") {
      val <- clustMixType::validation_kproto(
        method = "silhouette",
        data = dat$fullData,
        k = k_range,
        nstart = 2,
        verbose = FALSE
      )
      fit_kp <- clustMixType::kproto(dat$fullData, k = val$k_opt, nstart = 2, verbose = FALSE)
      t_elapsed <- (proc.time() - t_start)[["elapsed"]] * 1000
      ari <- calc_ari(dat$trueID, as.integer(fit_kp$cluster))
      data.frame(
        Method = m_name,
        Package = m_pkg,
        True_K = true_k,
        Predicted_K = val$k_opt,
        Criterion = "Silhouette Index",
        ARI = round(ari, 4),
        Time_ms = round(t_elapsed, 1),
        stringsAsFactors = FALSE
      )
    } else if (m == "varsellcm") {
      res_v <- VarSelLCM::VarSelCluster(
        x = dat$fullData, gvals = k_range, vbleSelec = FALSE, crit.varsel = "BIC", nbcores = 1
      )
      t_elapsed <- (proc.time() - t_start)[["elapsed"]] * 1000
      memb_v <- as.integer(VarSelLCM::fitted(res_v, type = "partition"))
      pred_k <- tryCatch({
        if (methods::.hasSlot(res_v, "model") && methods::.hasSlot(res_v@model, "g")) {
          as.integer(res_v@model@g)
        } else {
          length(unique(memb_v))
        }
      }, error = function(e) length(unique(memb_v)))
      ari <- calc_ari(dat$trueID, memb_v)
      data.frame(
        Method = m_name,
        Package = m_pkg,
        True_K = true_k,
        Predicted_K = pred_k,
        Criterion = "BIC / MICL",
        ARI = round(ari, 4),
        Time_ms = round(t_elapsed, 1),
        stringsAsFactors = FALSE
      )
    } else if (m == "flexmix_multinom") {
      con_cols <- colnames(dat$conVars)
      cat_cols <- colnames(dat$catVars)
      con_form <- stats::as.formula(paste0("cbind(", paste(con_cols, collapse = ", "), ") ~ 1"))
      con_mod <- flexmix::FLXMCmvnorm(con_form, diagonal = TRUE)
      cat_mods <- lapply(cat_cols, function(col) {
        flexmix::FLXMRmultinom(stats::as.formula(paste0(col, " ~ 1")))
      })
      m_step <- flexmix::stepFlexmix(
        stats::as.formula(paste0(con_cols[1], " ~ 1")),
        data = dat$fullData,
        k = k_range,
        nrep = 1,
        model = c(list(con_mod), cat_mods),
        control = list(iter.max = 20, minprior = 0.05, verbose = 0)
      )
      best_m <- flexmix::getModel(m_step, which = "BIC")
      t_elapsed <- (proc.time() - t_start)[["elapsed"]] * 1000
      ari <- calc_ari(dat$trueID, as.integer(flexmix::clusters(best_m)))
      data.frame(
        Method = m_name,
        Package = m_pkg,
        True_K = true_k,
        Predicted_K = best_m@k,
        Criterion = "BIC (Multinomial)",
        ARI = round(ari, 4),
        Time_ms = round(t_elapsed, 1),
        stringsAsFactors = FALSE
      )
    } else if (m == "flexmix_binary") {
      con_mat <- safe_scale(dat$conVars)
      cat_dummy <- stats::model.matrix(~ . - 1, data = dat$catVars)
      df_comb <- as.data.frame(cbind(con_mat, cat_dummy))
      con_cols <- colnames(con_mat)
      bin_cols <- colnames(cat_dummy)
      con_form <- stats::as.formula(paste0("cbind(", paste(con_cols, collapse = ", "), ") ~ 1"))
      bin_form <- stats::as.formula(paste0("cbind(", paste(bin_cols, collapse = ", "), ") ~ 1"))
      m_step <- flexmix::stepFlexmix(
        stats::as.formula(paste0(con_cols[1], " ~ 1")),
        data = df_comb,
        k = k_range,
        nrep = 1,
        model = list(
          flexmix::FLXMCmvnorm(con_form, diagonal = TRUE),
          flexmix::FLXMCmvbinary(bin_form)
        ),
        control = list(iter.max = 20, minprior = 0.05, verbose = 0)
      )
      best_m <- flexmix::getModel(m_step, which = "BIC")
      t_elapsed <- (proc.time() - t_start)[["elapsed"]] * 1000
      ari <- calc_ari(dat$trueID, as.integer(flexmix::clusters(best_m)))
      data.frame(
        Method = m_name,
        Package = m_pkg,
        True_K = true_k,
        Predicted_K = best_m@k,
        Criterion = "BIC (Binary Levels)",
        ARI = round(ari, 4),
        Time_ms = round(t_elapsed, 1),
        stringsAsFactors = FALSE
      )
    }
  }, error = function(e) {
    data.frame(
      Method = m_name,
      Package = m_pkg,
      True_K = true_k,
      Predicted_K = NA,
      Criterion = paste("Error:", substr(e$message, 1, 24)),
      ARI = NA,
      Time_ms = NA,
      stringsAsFactors = FALSE
    )
  })
}

run_parallel_jobs <- function(method_keys, runner_fn, use_parallel = TRUE, num_cores = 4,
                              cluster_obj = NULL, progress_fn = NULL, ...) {
  is_webr <- exists("webr", envir = .GlobalEnv) || Sys.getenv("WEBR") == "1"
  n_methods <- length(method_keys)
  if (!use_parallel || is_webr || num_cores <= 1 || n_methods <= 1) {
    t_start <- proc.time()
    res <- vector("list", n_methods)
    for (i in seq_along(method_keys)) {
      m <- method_keys[i]
      if (is.function(progress_fn)) {
        progress_fn(i, n_methods, method_meta[[m]]$name)
      }
      res[[i]] <- runner_fn(m, ...)
    }
    t_elapsed <- (proc.time() - t_start)[["elapsed"]] * 1000
    attr(res, "wall_clock_ms") <- t_elapsed
    return(res)
  }

  n_workers <- min(n_methods, as.integer(num_cores))
  cl <- if (!is.null(cluster_obj)) cluster_obj else tryCatch(parallel::makeCluster(n_workers), error = function(e) NULL)
  is_temp_cl <- is.null(cluster_obj) && !is.null(cl)

  if (is.null(cl)) {
    t_start <- proc.time()
    res <- vector("list", n_methods)
    for (i in seq_along(method_keys)) {
      m <- method_keys[i]
      if (is.function(progress_fn)) {
        progress_fn(i, n_methods, method_meta[[m]]$name)
      }
      res[[i]] <- runner_fn(m, ...)
    }
    t_elapsed <- (proc.time() - t_start)[["elapsed"]] * 1000
    attr(res, "wall_clock_ms") <- t_elapsed
    return(res)
  }

  if (is_temp_cl) {
    on.exit(try(parallel::stopCluster(cl), silent = TRUE), add = TRUE)
  }

  args_list <- list(...)
  worker_exec <- function(m, args) {
    res <- do.call(runner_fn, c(list(m = m), args))
    list(m = m, result = res)
  }

  parallel::clusterExport(
    cl,
    varlist = c(
      "has_pkg", "safe_scale", "calc_ari", "calc_misclass_error",
      "method_meta", "runner_fn", "args_list", "worker_exec"
    ),
    envir = environment()
  )

  t_start <- proc.time()
  submitted <- min(n_workers, n_methods)
  for (i in seq_len(submitted)) {
    parallel:::sendCall(
      cl[[i]],
      worker_exec,
      list(m = method_keys[i], args = args_list)
    )
  }

  res_map <- vector("list", n_methods)
  names(res_map) <- method_keys

  for (i in seq_len(n_methods)) {
    worker_res <- parallel:::recvOneResult(cl)
    m_done <- worker_res$value$m
    res_map[[m_done]] <- worker_res$value$result

    if (is.function(progress_fn)) {
      m_name <- if (m_done %in% names(method_meta)) method_meta[[m_done]]$name else m_done
      progress_fn(i, n_methods, sprintf("Completed %s", m_name))
    }

    if (submitted < n_methods) {
      submitted <- submitted + 1
      parallel:::sendCall(
        cl[[worker_res$node]],
        worker_exec,
        list(m = method_keys[submitted], args = args_list)
      )
    }
  }

  t_elapsed <- (proc.time() - t_start)[["elapsed"]] * 1000
  res <- unname(res_map[method_keys])
  attr(res, "wall_clock_ms") <- t_elapsed
  res
}

# ------------------------------------------------------------------------------
# UI Layout
# ------------------------------------------------------------------------------
ui <- fluidPage(
  theme = if (has_pkg("bslib")) bslib::bs_theme(version = 5, bootswatch = "flatly") else NULL,

  tags$head(
    tags$style(HTML("
      .shiny-notification {
        position: fixed;
        top: 20px;
        right: 20px;
        width: 340px;
        border-radius: 8px !important;
        box-shadow: 0 6px 16px rgba(0,0,0,0.2) !important;
        border-left: 6px solid #18bc9c !important;
        font-weight: 500 !important;
        z-index: 9999;
      }
    "))
  ),

  titlePanel(
    tags$div(
      tags$h2("Mixed-Type Clustering Horse-Race", style = "margin-bottom: 2px; font-weight: 700;"),
      tags$p(
        "Interactive benchmark comparing mixed-data clustering techniques.",
        style = "color: #6c757d; font-size: 1.05rem;"
      )
    ),
    windowTitle = "Mixed-Type Clustering Horse-Race"
  ),

  sidebarLayout(
    sidebarPanel(
      width = 4,
      tags$h4("Simulation Settings", style = "font-weight: 600;"),

      sliderInput(
        "log_n",
        "Number of Observations log10(N):",
        min = 3.0,
        max = 5.0,
        value = 3.0,
        step = 0.1
      ),
      tags$div(
        id = "n_obs_badge",
        class = "shiny-html-output",
        style = "margin-top: -10px; margin-bottom: 12px;",
        tags$div(
          tags$span(
            class = "badge bg-light text-dark border",
            style = "font-size: 0.88rem; padding: 4px 8px; font-weight: 600;",
            "Selected N = 1,000 (10^3.0)"
          ),
          tags$span(
            style = "margin-left: 6px; font-size: 0.80rem; color: #6c757d;",
            "(Logarithmic: 10^3.0 = 1,000 to 10^5.0 = 100,000)"
          )
        )
      ),
      fluidRow(
        column(6, sliderInput("p_con", "Num. Continuous Vars.", min = 2, max = 30, value = 10, step = 1)),
        column(6, sliderInput("p_cat", "Num. Categorical Vars.", min = 2, max = 30, value = 10, step = 1))
      ),
      fluidRow(
        column(6, sliderInput("k_clusters", "Num. True Clusters:", min = 2, max = 10, value = 4, step = 1)),
        column(6, sliderInput("num_levels", "Num. Levels per Cat. Var:", min = 2, max = 20, value = 10, step = 1))
      ),
      fluidRow(
        column(6, sliderInput("separation", "Avg. Cluster Separation", min = 0.5, max = 4.0, value = 2.0, step = 0.5)),
        column(6, numericInput("rand_seed", "Random Seed:", value = 42, min = 1))
      ),

      tags$hr(),
      tags$h4("Techniques to Compare", style = "font-weight: 600;"),
      checkboxGroupInput(
        "methods",
        label = NULL,
        choices = c(
          "KAMILA (kamila)" = "kamila",
          "Gower's Dist + PAM (cluster)" = "gower_pam",
          "K-Prototypes (clustMixType)" = "kproto",
          "VarSelLCM (VarSelLCM)" = "varsellcm",
          "FlexMix Multinomial (flexmix)" = "flexmix_multinom",
          "FlexMix Binary Levels (flexmix)" = "flexmix_binary"
        ),
        selected = c(
          "kamila", "gower_pam", "kproto", "varsellcm",
          "flexmix_multinom", "flexmix_binary"
        )
      ),

      tags$hr(),
      tags$h4("Compute & Parallelism", style = "font-weight: 600;"),
      checkboxInput(
        "use_parallel",
        "Enable Multi-Core Parallel Execution",
        value = TRUE
      ),
      conditionalPanel(
        condition = "input.use_parallel",
        sliderInput(
          "num_cores",
          "Worker Cores:",
          min = 1,
          max = detected_cores_count,
          value = min(4, detected_cores_count),
          step = 1
        )
      )
    ),

    mainPanel(
      width = 8,
      tabsetPanel(
        id = "main_tabs",

        tabPanel(
          "Fixed-K Benchmark",
          tags$br(),
          tags$div(
            class = "alert alert-info",
            "Evaluate ARI (Adjusted Rand Index), misclassification error, and runtime for a fixed number of clusters."
          ),
          actionButton(
            "btn_run_fixed",
            "Run Fixed-K Benchmark",
            class = "btn-primary btn-md",
            style = "margin-top: 4px; margin-bottom: 16px; font-weight: 600;"
          ),
          tags$h4("Performance Summary", style = "font-weight: 600; margin-top: 10px;"),
          uiOutput("fixed_timing_summary"),
          tableOutput("benchmark_table"),
          tags$hr(),
          tags$h4("Visual Performance Comparison", style = "font-weight: 600;"),
          plotOutput("benchmark_plot", height = "420px"),
          tags$hr(),
          tags$h4("Executable R Code Snippets", style = "font-weight: 600;"),
          selectInput(
            "code_method_fixed",
            "Select Technique for R Code:",
            choices = c(
              "KAMILA (kamila)" = "kamila",
              "Gower + PAM (cluster)" = "gower_pam",
              "K-Prototypes (clustMixType)" = "kproto",
              "VarSelLCM (VarSelLCM)" = "varsellcm",
              "FlexMix (Multinomial)" = "flexmix_multinom",
              "FlexMix (Binary Levels)" = "flexmix_binary"
            ),
            selected = "kamila"
          ),
          verbatimTextOutput("code_fixed_out")
        ),

        tabPanel(
          "Cluster Number Selection (K in 2:10)",
          tags$br(),
          tags$div(
            class = "alert alert-success",
            tags$strong("Model Selection: "),
            "Simulates unknown cluster count over candidate K in [2, 10]. Shows predicted cluster count and criterion."
          ),
          actionButton(
            "btn_run_select",
            "Run Cluster Selection (K in 2:10)",
            class = "btn-primary btn-md",
            style = "margin-top: 4px; margin-bottom: 16px; font-weight: 600;"
          ),
          tags$h4("Cluster Selection Results", style = "font-weight: 600; margin-top: 10px;"),
          uiOutput("selection_timing_summary"),
          tableOutput("selection_table"),
          tags$hr(),
          tags$h4("Selected K Comparison Plot", style = "font-weight: 600;"),
          plotOutput("selection_plot", height = "400px"),
          tags$hr(),
          tags$h4("Executable R Code Snippets (Model Selection)", style = "font-weight: 600;"),
          selectInput(
            "code_method_sel",
            "Select Technique for Model Selection Code:",
            choices = c(
              "KAMILA (kamila)" = "kamila",
              "Gower + PAM (cluster)" = "gower_pam",
              "K-Prototypes (clustMixType)" = "kproto",
              "VarSelLCM (VarSelLCM)" = "varsellcm",
              "FlexMix (Multinomial)" = "flexmix_multinom",
              "FlexMix (Binary Levels)" = "flexmix_binary"
            ),
            selected = "kamila"
          ),
          verbatimTextOutput("code_sel_out")
        ),

        tabPanel(
          "LDA Cluster Projection",
          tags$br(),
          tags$div(
            class = "alert alert-secondary",
            tags$strong("Linear Discriminant Analysis (LDA) Projection: "),
            "Displays the multi-dimensional mixed data projected onto the optimal discriminant subspace ",
            "defined by the true ground-truth clusters. Compare partitioning across all techniques."
          ),
          fluidRow(
            column(
              6,
              selectInput(
                "proj_view",
                "Projection Layout:",
                choices = c(
                  "Multi-Method Comparison Grid" = "grid",
                  "Single Method Focus" = "single"
                ),
                selected = "grid"
              )
            ),
            column(
              6,
              conditionalPanel(
                condition = "input.proj_view == 'single'",
                selectInput(
                  "proj_single_method",
                  "Select Technique to Display:",
                  choices = c(
                    "True Ground Truth" = "true",
                    "KAMILA (kamila)" = "kamila",
                    "Gower + PAM (cluster)" = "gower_pam",
                    "K-Prototypes (clustMixType)" = "kproto",
                    "VarSelLCM (VarSelLCM)" = "varsellcm",
                    "FlexMix (Multinomial)" = "flexmix_multinom",
                    "FlexMix (Binary Levels)" = "flexmix_binary"
                  ),
                  selected = "kamila"
                )
              )
            )
          ),
          plotOutput("lda_cluster_plot", height = "600px"),
          tags$hr(),
          tags$h4("Executable R Code Snippet (Clustering & LDA Projection)", style = "font-weight: 600;"),
          verbatimTextOutput("code_lda_out")
        ),

        tabPanel(
          "Environment & Hardware",
          tags$br(),
          tags$h4("Runtime Environment & Specs", style = "font-weight: 600;"),
          tableOutput("env_info_table"),
          tags$h4("Package Versions & Availability", style = "font-weight: 600; margin-top: 20px;"),
          tableOutput("pkg_version_table")
        )
      )
    )
  )
)

# ------------------------------------------------------------------------------
# Server Logic
# ------------------------------------------------------------------------------
server <- function(input, output, session) {

  # Session-Level Worker Cluster Lifecycle Management
  session_cluster <- reactiveVal(NULL)
  current_cluster_size <- reactiveVal(0L)
  fixed_wall_time <- reactiveVal(NULL)
  selection_wall_time <- reactiveVal(NULL)

  get_or_create_cluster <- function(n_workers) {
    cl <- session_cluster()
    cur_size <- current_cluster_size()
    if (!is.null(cl) && cur_size == n_workers) {
      return(cl)
    }
    if (!is.null(cl)) {
      try(parallel::stopCluster(cl), silent = TRUE)
      session_cluster(NULL)
      current_cluster_size(0L)
    }
    new_cl <- tryCatch(parallel::makeCluster(n_workers), error = function(e) NULL)
    if (is.null(new_cl)) return(NULL)

    parallel::clusterEvalQ(new_cl, {
      suppressPackageStartupMessages({
        library(kamila)
        library(cluster)
        if (requireNamespace("clustMixType", quietly = TRUE)) library(clustMixType)
        if (requireNamespace("VarSelLCM", quietly = TRUE)) library(VarSelLCM)
        if (requireNamespace("flexmix", quietly = TRUE)) library(flexmix)
        if (requireNamespace("mvtnorm", quietly = TRUE)) library(mvtnorm)
        if (requireNamespace("mclust", quietly = TRUE)) library(mclust)
      })
    })
    session_cluster(new_cl)
    current_cluster_size(n_workers)
    new_cl
  }

  session$onSessionEnded(function() {
    cl <- session_cluster()
    if (!is.null(cl)) {
      try(parallel::stopCluster(cl), silent = TRUE)
    }
  })

  output$n_obs_badge <- renderUI({
    val <- if (is.null(input$log_n) || !is.numeric(input$log_n)) 3.0 else input$log_n
    n_val <- as.integer(round(10^val))
    tags$div(
      tags$span(
        class = "badge bg-light text-dark border",
        style = "font-size: 0.88rem; padding: 4px 8px; font-weight: 600;",
        sprintf("Selected N = %s (10^%.1f)", format(n_val, big.mark = ","), val)
      ),
      tags$span(
        style = "margin-left: 6px; font-size: 0.80rem; color: #6c757d;",
        "(Logarithmic: 10^3.0 = 1,000 to 10^5.0 = 100,000)"
      )
    )
  })

  # Dataset Generator
  sim_data <- reactive({
    input$btn_run_fixed
    input$btn_run_select
    isolate({
      n_calc <- if (!is.null(input$log_n)) as.integer(round(10^input$log_n)) else 1000L
      generate_synthetic_mixed_data(
        n = n_calc,
        p_con = input$p_con,
        p_cat = input$p_cat,
        k = input$k_clusters,
        separation = input$separation,
        num_levels = input$num_levels,
        seed = input$rand_seed
      )
    })
  })

  # ----------------------------------------------------------------------------
  # Fixed-K Benchmark
  # ----------------------------------------------------------------------------
  benchmark_results <- eventReactive(list(input$btn_run_fixed, input$rand_seed), {
    dat <- sim_data()
    selected_methods <- input$methods
    k <- isolate(input$k_clusters)
    use_par <- isTRUE(input$use_parallel)
    n_cores <- if (!is.null(input$num_cores)) as.integer(input$num_cores) else 1

    valid_methods <- intersect(names(method_meta), selected_methods)
    if (length(valid_methods) == 0) {
      return(data.frame(Message = "No techniques selected"))
    }

    msg <- if (use_par && n_cores > 1) {
      sprintf(
        "Running Fixed-K Benchmark (Parallel across %d cores)...",
        min(n_cores, length(valid_methods))
      )
    } else {
      "Running Fixed-K Benchmark (Sequential)..."
    }

    worker_cl <- if (use_par && n_cores > 1) {
      get_or_create_cluster(min(n_cores, length(valid_methods)))
    } else {
      NULL
    }

    results_list <- withProgress(
      message = msg,
      detail = "Initializing benchmark...",
      value = 0.05,
      {
        progress_cb <- function(cur_idx, total_count, item_desc) {
          incProgress(
            amount = 0.9 / total_count,
            detail = sprintf("[%d/%d] %s", cur_idx, total_count, item_desc)
          )
        }
        run_parallel_jobs(
          method_keys = valid_methods,
          runner_fn = run_single_fixed_k,
          use_parallel = use_par,
          num_cores = n_cores,
          cluster_obj = worker_cl,
          progress_fn = progress_cb,
          dat = dat,
          k = k
        )
      }
    )

    fixed_wall_time(attr(results_list, "wall_clock_ms"))
    do.call(rbind, results_list)
  }, ignoreNULL = FALSE)

  # ----------------------------------------------------------------------------
  # Cluster Selection Benchmark (K in 2:10)
  # ----------------------------------------------------------------------------
  selection_results <- eventReactive(input$btn_run_select, {
    dat <- sim_data()
    selected_methods <- input$methods
    true_k <- isolate(input$k_clusters)
    use_par <- isTRUE(input$use_parallel)
    n_cores <- if (!is.null(input$num_cores)) as.integer(input$num_cores) else 1

    valid_methods <- intersect(names(method_meta), selected_methods)
    if (length(valid_methods) == 0) {
      return(data.frame(Message = "No techniques selected"))
    }

    msg <- if (use_par && n_cores > 1) {
      sprintf(
        "Running Cluster Selection (Parallel across %d cores)...",
        min(n_cores, length(valid_methods))
      )
    } else {
      "Running Cluster Selection (Sequential)..."
    }

    worker_cl <- if (use_par && n_cores > 1) {
      get_or_create_cluster(min(n_cores, length(valid_methods)))
    } else {
      NULL
    }

    results_list <- withProgress(
      message = msg,
      detail = "Evaluating candidate cluster counts...",
      value = 0.05,
      {
        progress_cb <- function(cur_idx, total_count, item_desc) {
          incProgress(
            amount = 0.9 / total_count,
            detail = sprintf("[%d/%d] %s", cur_idx, total_count, item_desc)
          )
        }
        ps_cores_arg <- if (use_par && length(valid_methods) == 1) n_cores else 1
        run_parallel_jobs(
          method_keys = valid_methods,
          runner_fn = run_single_selection,
          use_parallel = use_par,
          num_cores = n_cores,
          cluster_obj = worker_cl,
          progress_fn = progress_cb,
          dat = dat,
          true_k = true_k,
          ps_cores = ps_cores_arg
        )
      }
    )

    selection_wall_time(attr(results_list, "wall_clock_ms"))
    valid_rows <- results_list[!sapply(results_list, is.null)]
    if (length(valid_rows) == 0) {
      return(data.frame(Message = "No selection results returned"))
    }
    do.call(rbind, valid_rows)
  })

  # ----------------------------------------------------------------------------
  # Multi-Method Cluster Membership Cache for Projection
  # ----------------------------------------------------------------------------
  cluster_assignments <- reactive({
    dat <- sim_data()
    k <- isolate(input$k_clusters)
    selected_methods <- input$methods
    membs <- list(true = dat$trueID)

    # KAMILA
    if ("kamila" %in% selected_methods && has_pkg("kamila")) {
      kam_res <- tryCatch({
        as.integer(kamila::kamila(
          dat$conVars,
          dat$catVars,
          numClust = k,
          numInit = 3,
          calcNumClust = "none"
        )$finalMemb)
      }, error = function(e) NULL)
      if (!is.null(kam_res)) membs$kamila <- kam_res
    }

    # Gower + PAM
    if ("gower_pam" %in% selected_methods && has_pkg("cluster")) {
      pam_res <- tryCatch({
        g_dist <- cluster::daisy(dat$fullData, metric = "gower")
        as.integer(cluster::pam(g_dist, k = k, diss = TRUE)$clustering)
      }, error = function(e) NULL)
      if (!is.null(pam_res)) membs$gower_pam <- pam_res
    }

    # K-Prototypes
    if ("kproto" %in% selected_methods && has_pkg("clustMixType")) {
      kp_res <- tryCatch({
        as.integer(clustMixType::kproto(dat$fullData, k = k, nstart = 2, verbose = FALSE)$cluster)
      }, error = function(e) NULL)
      if (!is.null(kp_res)) membs$kproto <- kp_res
    }

    # VarSelLCM
    if ("varsellcm" %in% selected_methods && has_pkg("VarSelLCM")) {
      v_res <- tryCatch({
        v_fit <- VarSelLCM::VarSelCluster(
          x = dat$fullData, gvals = k, vbleSelec = FALSE, crit.varsel = "BIC", nbcores = 1
        )
        as.integer(VarSelLCM::fitted(v_fit, type = "partition"))
      }, error = function(e) NULL)
      if (!is.null(v_res)) membs$varsellcm <- v_res
    }

    # FlexMix Multinomial
    if ("flexmix_multinom" %in% selected_methods && has_pkg("flexmix")) {
      f_res <- tryCatch({
        con_cols <- colnames(dat$conVars)
        cat_cols <- colnames(dat$catVars)
        con_form <- stats::as.formula(paste0("cbind(", paste(con_cols, collapse = ", "), ") ~ 1"))
        con_mod <- flexmix::FLXMCmvnorm(con_form, diagonal = TRUE)
        cat_mods <- lapply(cat_cols, function(col) {
          flexmix::FLXMRmultinom(stats::as.formula(paste0(col, " ~ 1")))
        })
        m <- flexmix::flexmix(
          stats::as.formula(paste0(con_cols[1], " ~ 1")),
          data = dat$fullData,
          k = k,
          model = c(list(con_mod), cat_mods),
          control = list(iter.max = 20, minprior = 0.05, verbose = 0)
        )
        as.integer(flexmix::clusters(m))
      }, error = function(e) NULL)
      if (!is.null(f_res)) membs$flexmix_multinom <- f_res
    }

    # FlexMix Binary Levels
    if ("flexmix_binary" %in% selected_methods && has_pkg("flexmix")) {
      fb_res <- tryCatch({
        con_mat <- safe_scale(dat$conVars)
        cat_dummy <- stats::model.matrix(~ . - 1, data = dat$catVars)
        df_comb <- as.data.frame(cbind(con_mat, cat_dummy))
        con_cols <- colnames(con_mat)
        bin_cols <- colnames(cat_dummy)
        con_form <- stats::as.formula(paste0("cbind(", paste(con_cols, collapse = ", "), ") ~ 1"))
        bin_form <- stats::as.formula(paste0("cbind(", paste(bin_cols, collapse = ", "), ") ~ 1"))
        m <- flexmix::flexmix(
          stats::as.formula(paste0(con_cols[1], " ~ 1")),
          data = df_comb,
          k = k,
          model = list(
            flexmix::FLXMCmvnorm(con_form, diagonal = TRUE),
            flexmix::FLXMCmvbinary(bin_form)
          ),
          control = list(iter.max = 20, minprior = 0.05, verbose = 0)
        )
        as.integer(flexmix::clusters(m))
      }, error = function(e) NULL)
      if (!is.null(fb_res)) membs$flexmix_binary <- fb_res
    }

    membs
  })

  # Method code to display name mapping
  method_name_map <- c(
    kamila = "KAMILA",
    gower_pam = "Gower + PAM",
    kproto = "K-Prototypes",
    varsellcm = "VarSelLCM",
    flexmix_multinom = "FlexMix (Multinomial)",
    flexmix_binary = "FlexMix (Binary Levels)"
  )

  # Render Tables & Plots
  output$fixed_timing_summary <- renderUI({
    wall_ms <- fixed_wall_time()
    res <- benchmark_results()
    if (is.null(wall_ms) || is.null(res) || !"Time_ms" %in% names(res)) return(NULL)
    active_methods <- method_name_map[input$methods]
    valid_res <- res[!is.na(res$Time_ms) & res$Method %in% active_methods, , drop = FALSE]
    if (nrow(valid_res) == 0) return(NULL)
    seq_sum <- sum(valid_res$Time_ms, na.rm = TRUE)
    speedup <- if (wall_ms > 0) round(seq_sum / wall_ms, 2) else 1.0
    tags$div(
      style = "margin-bottom: 12px; font-size: 0.95rem; color: #495057;",
      tags$span(tags$strong("Total Wall-Clock Time: "), sprintf("%.1f ms", wall_ms)),
      tags$span(" | "),
      tags$span(tags$strong("Sequential Sum of Runtimes: "), sprintf("%.1f ms", seq_sum)),
      if (speedup > 1.05) {
        tags$span(
          class = "badge bg-success",
          style = "margin-left: 8px; font-size: 0.85rem; padding: 4px 8px;",
          sprintf("%.2fx Parallel Speedup", speedup)
        )
      } else {
        NULL
      }
    )
  })

  output$benchmark_table <- renderTable({
    res <- benchmark_results()
    if (is.null(res) || !"Method" %in% names(res) || nrow(res) == 0) return(res)
    active_methods <- method_name_map[input$methods]
    filtered_res <- res[res$Method %in% active_methods, , drop = FALSE]
    if (nrow(filtered_res) == 0) {
      return(data.frame(Message = "No selected methods to display"))
    }
    filtered_res
  }, striped = TRUE, hover = TRUE, bordered = TRUE)

  output$benchmark_plot <- renderPlot({
    res <- benchmark_results()
    if (is.null(res) || !"ARI" %in% names(res) || nrow(res) == 0) return(NULL)
    active_methods <- method_name_map[input$methods]
    valid_res <- res[!is.na(res$ARI) & res$Method %in% active_methods, , drop = FALSE]
    if (nrow(valid_res) == 0) return(NULL)

    par(mfrow = c(1, 2), mar = c(10.5, 4.5, 3, 1))

    # ARI Plot
    barplot(
      valid_res$ARI,
      names.arg = valid_res$Method,
      col = "#2c3e50",
      main = "Adjusted Rand Index (Higher = Better)",
      ylab = "ARI Score",
      ylim = c(0, 1),
      las = 2,
      cex.names = 0.82
    )
    abline(h = seq(0, 1, 0.2), col = "gray80", lty = 2)

    # Timing Plot
    barplot(
      valid_res$Time_ms,
      names.arg = valid_res$Method,
      col = "#18bc9c",
      main = "Execution Time (Lower = Faster)",
      ylab = "Time (ms)",
      las = 2,
      cex.names = 0.82
    )
    abline(h = axTicks(2), col = "gray80", lty = 2)
  })

  output$selection_timing_summary <- renderUI({
    wall_ms <- selection_wall_time()
    res <- selection_results()
    if (is.null(wall_ms) || is.null(res) || !"Time_ms" %in% names(res)) return(NULL)
    active_methods <- method_name_map[input$methods]
    valid_res <- res[!is.na(res$Time_ms) & res$Method %in% active_methods, , drop = FALSE]
    if (nrow(valid_res) == 0) return(NULL)
    seq_sum <- sum(valid_res$Time_ms, na.rm = TRUE)
    speedup <- if (wall_ms > 0) round(seq_sum / wall_ms, 2) else 1.0
    tags$div(
      style = "margin-bottom: 12px; font-size: 0.95rem; color: #495057;",
      tags$span(tags$strong("Total Wall-Clock Time: "), sprintf("%.1f ms", wall_ms)),
      tags$span(" | "),
      tags$span(tags$strong("Sequential Sum of Runtimes: "), sprintf("%.1f ms", seq_sum)),
      if (speedup > 1.05) {
        tags$span(
          class = "badge bg-success",
          style = "margin-left: 8px; font-size: 0.85rem; padding: 4px 8px;",
          sprintf("%.2fx Parallel Speedup", speedup)
        )
      } else {
        NULL
      }
    )
  })

  output$selection_table <- renderTable({
    res <- selection_results()
    if (is.null(res) || !"Method" %in% names(res) || nrow(res) == 0) return(res)
    active_methods <- method_name_map[input$methods]
    filtered_res <- res[res$Method %in% active_methods, , drop = FALSE]
    if (nrow(filtered_res) == 0) {
      return(data.frame(Message = "No selected methods to display"))
    }
    filtered_res
  }, striped = TRUE, hover = TRUE, bordered = TRUE)

  output$selection_plot <- renderPlot({
    res <- selection_results()
    if (is.null(res) || !"Predicted_K" %in% names(res) || nrow(res) == 0) return(NULL)
    active_methods <- method_name_map[input$methods]
    valid_res <- res[!is.na(res$Predicted_K) & res$Method %in% active_methods, , drop = FALSE]
    if (nrow(valid_res) == 0) return(NULL)

    par(mar = c(10.5, 4.5, 3, 1))
    barplot(
      valid_res$Predicted_K,
      names.arg = valid_res$Method,
      col = "#3498db",
      main = "Predicted Number of Clusters (True K indicated by dashed line)",
      ylab = "Predicted K",
      ylim = c(0, 11),
      las = 2,
      cex.names = 0.85
    )
    abline(h = isolate(input$k_clusters), col = "red", lty = 2, lwd = 2)
    legend("topright", legend = paste("True K =", isolate(input$k_clusters)), col = "red", lty = 2, lwd = 2)
  })

  # ----------------------------------------------------------------------------
  # LDA Cluster Projection Plot (Subsampled 1000 Points with Multi-Method Colors)
  # ----------------------------------------------------------------------------
  output$lda_cluster_plot <- renderPlot({
    dat <- sim_data()
    all_membs <- cluster_assignments()
    k_true <- isolate(input$k_clusters)
    n_total <- nrow(dat$fullData)

    # Subsample 1,000 points (common across all panels and views)
    sample_size <- min(1000L, n_total)
    sub_idx <- if (n_total > sample_size) {
      set.seed(42)
      sample.int(n_total, size = sample_size)
    } else {
      seq_len(n_total)
    }

    # Construct design matrix for LDA projection using the subsample
    cat_mat <- model.matrix(~ ., data = dat$catVars[sub_idx, , drop = FALSE])[, -1, drop = FALSE]
    comb_mat <- cbind(scale(as.matrix(dat$conVars[sub_idx, , drop = FALSE])), scale(cat_mat))
    zv <- apply(comb_mat, 2, function(x) var(x, na.rm = TRUE) == 0 || is.na(var(x)))
    comb_mat <- comb_mat[, !zv, drop = FALSE]

    df_lda <- as.data.frame(comb_mat)
    df_lda$class <- as.factor(dat$trueID[sub_idx])

    # Compute LDA using true ground truth classes on subsample
    lda_fit <- MASS::lda(class ~ ., data = df_lda, tol = 1e-4)
    lda_pred <- predict(lda_fit, newdata = df_lda)
    lda_x <- lda_pred$x

    # Coordinate setup: If K == 2, LD1 + 1st Con variable; if K >= 3, LD1 + LD2
    if (ncol(lda_x) >= 2) {
      coords_x <- lda_x[, 1]
      coords_y <- lda_x[, 2]
      xlab_txt <- "Linear Discriminant 1 (LD1)"
      ylab_txt <- "Linear Discriminant 2 (LD2)"
    } else {
      coords_x <- lda_x[, 1]
      coords_y <- dat$conVars[sub_idx, 1]
      xlab_txt <- "Linear Discriminant 1 (LD1)"
      ylab_txt <- "Continuous Feature 1 (Con_1)"
    }

    palette <- c(
      "#e74c3c", "#3498db", "#2ecc71", "#f39c12", "#9b59b6",
      "#1abc9c", "#e67e22", "#34495e", "#d35400", "#16a085"
    )

    method_labels <- c(
      true = "True Ground Truth",
      kamila = "KAMILA (kamila)",
      gower_pam = "Gower + PAM (cluster)",
      kproto = "K-Prototypes (clustMixType)",
      varsellcm = "VarSelLCM",
      flexmix_multinom = "FlexMix (Multinomial)",
      flexmix_binary = "FlexMix (Binary Levels)"
    )

    sub_note <- if (n_total > sample_size) sprintf(" (Subsample N = %d of %d)", sample_size, n_total) else ""

    if (input$proj_view == "single") {
      # Single Focused Plot
      target_m <- input$proj_single_method
      full_memb <- if (target_m %in% names(all_membs)) all_membs[[target_m]] else dat$trueID
      cur_memb <- full_memb[sub_idx]
      title_str <- method_labels[target_m]
      if (is.na(title_str)) title_str <- target_m

      par(mar = c(5, 5, 4, 2))
      plot(
        coords_x, coords_y,
        col = palette[((cur_memb - 1) %% length(palette)) + 1],
        pch = 19,
        cex = 1.0,
        main = paste0("LDA Discriminant Projection: ", title_str, sub_note),
        xlab = xlab_txt,
        ylab = ylab_txt
      )
      grid()
      legend(
        "topright",
        legend = paste("Cluster", seq_len(k_true)),
        col = palette[((seq_len(k_true) - 1) %% length(palette)) + 1],
        pch = 19,
        bg = "white"
      )
    } else {
      # Multi-Method Grid View
      avail_keys <- names(all_membs)
      n_plots <- length(avail_keys)
      n_cols <- if (n_plots <= 2) 2 else if (n_plots <= 4) 2 else 3
      n_rows <- ceiling(n_plots / n_cols)

      par(mfrow = c(n_rows, n_cols), mar = c(4, 4, 3, 1))

      for (m_key in avail_keys) {
        full_memb <- all_membs[[m_key]]
        cur_memb <- full_memb[sub_idx]
        m_title <- method_labels[m_key]
        if (is.na(m_title)) m_title <- m_key

        plot(
          coords_x, coords_y,
          col = palette[((cur_memb - 1) %% length(palette)) + 1],
          pch = if (m_key == "true") 19 else 17,
          cex = 0.8,
          main = paste0(m_title, sub_note),
          xlab = xlab_txt,
          ylab = ylab_txt
        )
        grid()
      }
    }
  })

  # ----------------------------------------------------------------------------
  # Executable Code Snippets
  # ----------------------------------------------------------------------------
  output$code_fixed_out <- renderText({
    m <- input$code_method_fixed
    k <- isolate(input$k_clusters)
    if (is.null(m) || m == "kamila") {
      sprintf(
        paste0(
          "# --- KAMILA Clustering ---\n",
          "library(kamila)\n\n",
          "# dat$conVars: continuous variables (data.frame)\n",
          "# dat$catVars: categorical variables (data.frame of factors)\n",
          "fit <- kamila::kamila(\n",
          "  conVar = dat$conVars,\n",
          "  catFactor = dat$catVars,\n",
          "  numClust = %d,\n",
          "  numInit = 5,\n",
          "  maxIter = 25,\n",
          "  calcNumClust = \"none\"\n",
          ")\n",
          "clusters <- fit$finalMemb\n"
        ),
        k
      )
    } else if (m == "gower_pam") {
      sprintf(
        paste0(
          "# --- Gower Distance + PAM (Partitioning Around Medoids) ---\n",
          "library(cluster)\n\n",
          "# dat$fullData: mixed-type dataset\n",
          "g_dist <- cluster::daisy(dat$fullData, metric = \"gower\")\n",
          "fit <- cluster::pam(g_dist, k = %d, diss = TRUE)\n",
          "clusters <- fit$clustering\n"
        ),
        k
      )
    } else if (m == "kproto") {
      sprintf(
        paste0(
          "# --- K-Prototypes for Mixed-Type Data ---\n",
          "library(clustMixType)\n\n",
          "fit <- clustMixType::kproto(\n",
          "  x = dat$fullData,\n",
          "  k = %d,\n",
          "  nstart = 3,\n",
          "  verbose = FALSE\n",
          ")\n",
          "clusters <- fit$cluster\n"
        ),
        k
      )
    } else if (m == "varsellcm") {
      sprintf(
        paste0(
          "# --- VarSelLCM (Latent Class Model with Variable Selection) ---\n",
          "library(VarSelLCM)\n\n",
          "fit <- VarSelLCM::VarSelCluster(\n",
          "  x = dat$fullData,\n",
          "  gvals = %d,\n",
          "  vbleSelec = FALSE,\n",
          "  crit.varsel = \"BIC\",\n",
          "  nbcores = 1\n",
          ")\n",
          "clusters <- VarSelLCM::fitted(fit, type = \"partition\")\n"
        ),
        k
      )
    } else if (m == "flexmix_multinom") {
      sprintf(
        paste0(
          "# --- FlexMix (Joint Gaussian & Multinomial Model) ---\n",
          "library(flexmix)\n\n",
          "con_cols <- colnames(dat$conVars)\n",
          "cat_cols <- colnames(dat$catVars)\n",
          "con_form <- as.formula(paste0(\"cbind(\", paste(con_cols, collapse = \", \"), \") ~ 1\"))\n",
          "con_mod <- flexmix::FLXMCmvnorm(con_form, diagonal = TRUE)\n",
          "cat_mods <- lapply(cat_cols, function(col) {\n",
          "  flexmix::FLXMRmultinom(as.formula(paste0(col, \" ~ 1\")))\n",
          "})\n",
          "fit <- flexmix::flexmix(\n",
          "  as.formula(paste0(con_cols[1], \" ~ 1\")),\n",
          "  data = dat$fullData,\n",
          "  k = %d,\n",
          "  model = c(list(con_mod), cat_mods),\n",
          "  control = list(iter.max = 25, minprior = 0.05, verbose = 0)\n",
          ")\n",
          "clusters <- flexmix::clusters(fit)\n"
        ),
        k
      )
    } else if (m == "flexmix_binary") {
      sprintf(
        paste0(
          "# --- FlexMix (Joint Gaussian & Binary Level Indicators) ---\n",
          "library(flexmix)\n\n",
          "con_mat <- scale(as.matrix(dat$conVars))\n",
          "cat_dummy <- model.matrix(~ . - 1, data = dat$catVars)\n",
          "df_comb <- as.data.frame(cbind(con_mat, cat_dummy))\n",
          "con_cols <- colnames(con_mat)\n",
          "bin_cols <- colnames(cat_dummy)\n",
          "con_form <- as.formula(paste0(\"cbind(\", paste(con_cols, collapse = \", \"), \") ~ 1\"))\n",
          "bin_form <- as.formula(paste0(\"cbind(\", paste(bin_cols, collapse = \", \"), \") ~ 1\"))\n",
          "fit <- flexmix::flexmix(\n",
          "  as.formula(paste0(con_cols[1], \" ~ 1\")),\n",
          "  data = df_comb,\n",
          "  k = %d,\n",
          "  model = list(\n",
          "    flexmix::FLXMCmvnorm(con_form, diagonal = TRUE),\n",
          "    flexmix::FLXMCmvbinary(bin_form)\n",
          "  ),\n",
          "  control = list(iter.max = 25, minprior = 0.05, verbose = 0)\n",
          ")\n",
          "clusters <- flexmix::clusters(fit)\n"
        ),
        k
      )
    }
  })

  output$code_sel_out <- renderText({
    m <- input$code_method_sel
    k <- isolate(input$k_clusters)
    k_max <- min(10, max(5, k + 2))
    k_range_str <- sprintf("2:%d", k_max)

    if (is.null(m) || m == "kamila") {
      cores_snippet <- if (isTRUE(input$use_parallel) && !is.null(input$num_cores) && input$num_cores > 1) {
        sprintf(",\n  numCores = %d", as.integer(input$num_cores))
      } else {
        ""
      }
      sprintf(
        paste0(
          "# --- KAMILA Cluster Count Selection (Prediction Strength) ---\n",
          "library(kamila)\n\n",
          "fit <- kamila::kamila(\n",
          "  conVar = dat$conVars,\n",
          "  catFactor = dat$catVars,\n",
          "  numClust = %s,\n",
          "  numInit = 3,\n",
          "  calcNumClust = \"ps\",\n",
          "  numPredStrCvRun = 5,\n",
          "  predStrThresh = 0.6%s\n",
          ")\n",
          "best_k <- fit$nClust$bestNClust\n"
        ),
        k_range_str,
        cores_snippet
      )
    } else if (m == "gower_pam") {
      sprintf(
        paste0(
          "# --- Gower + PAM (Average Silhouette Width Selection) ---\n",
          "library(cluster)\n\n",
          "g_dist <- cluster::daisy(dat$fullData, metric = \"gower\")\n",
          "sils <- sapply(%s, function(ki) {\n",
          "  cluster::pam(g_dist, k = ki, diss = TRUE)$silinfo$avg.width\n",
          "})\n",
          "best_k <- (%s)[which.max(sils)]\n"
        ),
        k_range_str,
        k_range_str
      )
    } else if (m == "kproto") {
      sprintf(
        paste0(
          "# --- K-Prototypes Validation Index (Silhouette Selection) ---\n",
          "library(clustMixType)\n\n",
          "val <- clustMixType::validation_kproto(\n",
          "  method = \"silhouette\",\n",
          "  data = dat$fullData,\n",
          "  k = %s,\n",
          "  nstart = 2,\n",
          "  verbose = FALSE\n",
          ")\n",
          "best_k <- val$k_opt\n"
        ),
        k_range_str
      )
    } else if (m == "varsellcm") {
      sprintf(
        paste0(
          "# --- VarSelLCM (BIC / MICL Information Criterion Selection) ---\n",
          "library(VarSelLCM)\n\n",
          "fit <- VarSelLCM::VarSelCluster(\n",
          "  x = dat$fullData,\n",
          "  gvals = %s,\n",
          "  vbleSelec = FALSE,\n",
          "  crit.varsel = \"BIC\",\n",
          "  nbcores = 1\n",
          ")\n",
          "best_k <- fit@model@g\n"
        ),
        k_range_str
      )
    } else if (m == "flexmix_multinom") {
      sprintf(
        paste0(
          "# --- FlexMix (stepFlexmix BIC Multinomial Selection) ---\n",
          "library(flexmix)\n\n",
          "con_cols <- colnames(dat$conVars)\n",
          "cat_cols <- colnames(dat$catVars)\n",
          "con_form <- as.formula(paste0(\"cbind(\", paste(con_cols, collapse = \", \"), \") ~ 1\"))\n",
          "con_mod <- flexmix::FLXMCmvnorm(con_form, diagonal = TRUE)\n",
          "cat_mods <- lapply(cat_cols, function(col) {\n",
          "  flexmix::FLXMRmultinom(as.formula(paste0(col, \" ~ 1\")))\n",
          "})\n",
          "m_step <- flexmix::stepFlexmix(\n",
          "  as.formula(paste0(con_cols[1], \" ~ 1\")),\n",
          "  data = dat$fullData,\n",
          "  k = %s,\n",
          "  nrep = 1,\n",
          "  model = c(list(con_mod), cat_mods),\n",
          "  control = list(iter.max = 20, minprior = 0.05, verbose = 0)\n",
          ")\n",
          "best_model <- flexmix::getModel(m_step, which = \"BIC\")\n",
          "best_k <- best_model@k\n"
        ),
        k_range_str
      )
    } else if (m == "flexmix_binary") {
      sprintf(
        paste0(
          "# --- FlexMix (stepFlexmix BIC Binary Levels Selection) ---\n",
          "library(flexmix)\n\n",
          "con_mat <- scale(as.matrix(dat$conVars))\n",
          "cat_dummy <- model.matrix(~ . - 1, data = dat$catVars)\n",
          "df_comb <- as.data.frame(cbind(con_mat, cat_dummy))\n",
          "con_cols <- colnames(con_mat)\n",
          "bin_cols <- colnames(cat_dummy)\n",
          "con_form <- as.formula(paste0(\"cbind(\", paste(con_cols, collapse = \", \"), \") ~ 1\"))\n",
          "bin_form <- as.formula(paste0(\"cbind(\", paste(bin_cols, collapse = \", \"), \") ~ 1\"))\n",
          "m_step <- flexmix::stepFlexmix(\n",
          "  as.formula(paste0(con_cols[1], \" ~ 1\")),\n",
          "  data = df_comb,\n",
          "  k = %s,\n",
          "  nrep = 1,\n",
          "  model = list(\n",
          "    flexmix::FLXMCmvnorm(con_form, diagonal = TRUE),\n",
          "    flexmix::FLXMCmvbinary(bin_form)\n",
          "  ),\n",
          "  control = list(iter.max = 20, minprior = 0.05, verbose = 0)\n",
          ")\n",
          "best_model <- flexmix::getModel(m_step, which = \"BIC\")\n",
          "best_k <- best_model@k\n"
        ),
        k_range_str
      )
    }
  })

  output$code_lda_out <- renderText({
    paste0(
      "# --- 1. Compute Ground-Truth Linear Discriminant Analysis (LDA) ---\n",
      "library(MASS)\n\n",
      "# Subsample 1,000 points if N > 1,000\n",
      "n_total <- nrow(dat$fullData)\n",
      "sub_idx <- if (n_total > 1000) sample.int(n_total, 1000) else seq_len(n_total)\n",
      "cat_mat <- model.matrix(~ ., data = dat$catVars[sub_idx, , drop = FALSE])[, -1, drop = FALSE]\n",
      "comb_mat <- cbind(scale(as.matrix(dat$conVars[sub_idx, , drop = FALSE])), scale(cat_mat))\n",
      "df_lda <- as.data.frame(comb_mat)\n",
      "df_lda$class <- as.factor(dat$trueID[sub_idx])\n\n",
      "# Fit LDA model on true cluster classes\n",
      "lda_fit <- MASS::lda(class ~ ., data = df_lda, tol = 1e-4)\n",
      "lda_coords <- predict(lda_fit, newdata = df_lda)$x\n\n",
      "# --- 2. Project Points Colored by Predicted Cluster Partitions ---\n",
      "palette <- c(\"#e74c3c\", \"#3498db\", \"#2ecc71\", \"#f39c12\",\n",
      "             \"#9b59b6\", \"#1abc9c\", \"#e67e22\", \"#34495e\")\n",
      "plot(\n",
      "  lda_coords[, 1], lda_coords[, 2],\n",
      "  col = palette[clusters[sub_idx]],\n",
      "  pch = 19,\n",
      "  xlab = \"Linear Discriminant 1 (LD1)\",\n",
      "  ylab = \"Linear Discriminant 2 (LD2)\",\n",
      "  main = \"LDA Discriminant Subspace Projection (Subsample N = 1000)\"\n",
      ")\n",
      "grid()\n"
    )
  })

  # Environment Information
  output$env_info_table <- renderTable({
    is_webr <- exists("webr", envir = .GlobalEnv) || Sys.getenv("WEBR") == "1"
    data.frame(
      Property = c(
        "R Version",
        "Platform / Architecture",
        "Client Hardware Cores",
        "Execution Engine",
        "Hosting Architecture"
      ),
      Value = c(
        R.version.string,
        R.version$platform,
        as.character(parallel::detectCores(logical = FALSE)),
        if (is_webr) "WebR / WebAssembly (Client-Side in Browser)" else "Native R Runtime",
        if (is_webr) "Shinylive (Zero Server / Static GitHub Pages)" else "Standard Shiny Session"
      ),
      stringsAsFactors = FALSE
    )
  }, striped = TRUE, bordered = TRUE)

  # Package Versions
  output$pkg_version_table <- renderTable({
    pkgs <- c("kamila", "cluster", "clustMixType", "VarSelLCM", "flexmix", "mixtools", "mclust", "shiny")
    versions <- sapply(pkgs, function(p) {
      if (has_pkg(p)) as.character(packageVersion(p)) else "Not Installed (Optional)"
    })
    data.frame(
      Package = pkgs,
      Version = versions,
      stringsAsFactors = FALSE
    )
  }, striped = TRUE, bordered = TRUE)
}

# ------------------------------------------------------------------------------
# Launch Application
# ------------------------------------------------------------------------------
if (requireNamespace("shiny", quietly = TRUE)) {
  shinyApp(ui = ui, server = server)
}
