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
# Synthetic Mixed Data Generator (1,000 to 5,000 observations)
# ------------------------------------------------------------------------------
generate_synthetic_mixed_data <- function(n = 1000, p_con = 5, p_cat = 5, k = 3,
                                          separation = 2.0, num_levels = 4, seed = NULL) {
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
      .horse-status-badge {
        display: inline-block;
        padding: 4px 10px;
        border-radius: 12px;
        font-size: 0.85rem;
        font-weight: 600;
        margin: 2px;
      }
      .horse-running {
        animation: horsePulse 1.1s infinite ease-in-out;
      }
      @keyframes horsePulse {
        0% { transform: scale(1); opacity: 0.85; }
        50% { transform: scale(1.05); opacity: 1; }
        100% { transform: scale(1); opacity: 0.85; }
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
        "n_obs",
        "Number of Observations (N):",
        min = 1000,
        max = 5000,
        value = 1000,
        step = 250
      ),
      fluidRow(
        column(6, sliderInput("p_con", "Num. Continuous Vars.", min = 2, max = 15, value = 5, step = 1)),
        column(6, sliderInput("p_cat", "Num. Categorial Vars.", min = 2, max = 15, value = 5, step = 1))
      ),
      fluidRow(
        column(6, sliderInput("k_clusters", "Num. True Clusters:", min = 2, max = 5, value = 3, step = 1)),
        column(6, sliderInput("separation", "Avg. Cluster Separation", min = 0.5, max = 4.0, value = 2.0, step = 0.5))
      ),
      numericInput("rand_seed", "Random Seed:", value = 42, min = 1),

      tags$hr(),
      tags$h4("Techniques to Compare", style = "font-weight: 600;"),
      checkboxGroupInput(
        "methods",
        label = NULL,
        choices = c(
          "KAMILA (kamila)" = "kamila",
          "Gower's Dist + PAM (cluster)" = "gower_pam",
          "K-Prototypes (clustMixType)" = "kproto",
          "ClusPCAMix (clustrd)" = "cluspcamix",
          "VarSelLCM (VarSelLCM)" = "varsellcm",
          "FlexMix / Latent Class (flexmix)" = "flexmix"
        ),
        selected = c("kamila", "gower_pam", "kproto", "cluspcamix", "varsellcm", "flexmix")
      ),

      tags$hr(),
      uiOutput("live_status_ui")
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
          tableOutput("benchmark_table"),
          tags$hr(),
          tags$h4("Visual Performance Comparison", style = "font-weight: 600;"),
          plotOutput("benchmark_plot", height = "380px")
        ),

        tabPanel(
          "Cluster Number Selection (K in 2:5)",
          tags$br(),
          tags$div(
            class = "alert alert-success",
            tags$strong("Model Selection: "),
            "Simulates unknown cluster count over K in [2, 5]. Shows predicted cluster count and selection criterion."
          ),
          actionButton(
            "btn_run_select",
            "Run Cluster Selection (K in 2:5)",
            class = "btn-primary btn-md",
            style = "margin-top: 4px; margin-bottom: 16px; font-weight: 600;"
          ),
          tags$h4("Cluster Selection Results", style = "font-weight: 600; margin-top: 10px;"),
          tableOutput("selection_table"),
          tags$hr(),
          tags$h4("Selected K Comparison Plot", style = "font-weight: 600;"),
          plotOutput("selection_plot", height = "360px")
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
                    "ClusPCAMix (clustrd)" = "cluspcamix",
                    "VarSelLCM (VarSelLCM)" = "varsellcm",
                    "FlexMix (flexmix)" = "flexmix"
                  ),
                  selected = "kamila"
                )
              )
            )
          ),
          plotOutput("lda_cluster_plot", height = "600px")
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

  # Reactive execution progress state
  status_tracker <- reactiveVal(list())

  # Dataset Generator
  sim_data <- reactive({
    input$btn_run_fixed
    input$btn_run_select
    isolate({
      generate_synthetic_mixed_data(
        n = input$n_obs,
        p_con = input$p_con,
        p_cat = input$p_cat,
        k = input$k_clusters,
        separation = input$separation,
        seed = input$rand_seed
      )
    })
  })

  # Live progress badge panel
  output$live_status_ui <- renderUI({
    st <- status_tracker()
    if (length(st) == 0) return(NULL)

    tagList(
      tags$div(
        style = "padding: 8px; background: #f8f9fa; border-radius: 6px; border: 1px solid #dee2e6;",
        tags$h6("Execution Tracker", style = "font-weight: 700; margin-bottom: 6px;"),
        lapply(names(st), function(m) {
          info <- st[[m]]
          status_class <- if (info$state == "running") {
            "badge bg-warning text-dark horse-running"
          } else if (info$state == "done") {
            "badge bg-success"
          } else if (info$state == "error") {
            "badge bg-danger"
          } else {
            "badge bg-secondary"
          }
          tags$div(
            style = "margin-bottom: 4px; font-size: 0.85rem;",
            tags$span(class = status_class, info$label),
            tags$span(style = "margin-left: 6px; color: #495057;", info$detail)
          )
        })
      )
    )
  })

  # ----------------------------------------------------------------------------
  # Fixed-K Benchmark
  # ----------------------------------------------------------------------------
  benchmark_results <- eventReactive(list(input$btn_run_fixed, input$rand_seed), {
    dat <- sim_data()
    selected_methods <- input$methods
    k <- isolate(input$k_clusters)
    results <- list()

    method_meta <- list(
      kamila = list(name = "KAMILA", pkg = "kamila"),
      gower_pam = list(name = "Gower + PAM", pkg = "cluster"),
      kproto = list(name = "K-Prototypes", pkg = "clustMixType"),
      cluspcamix = list(name = "ClusPCAMix", pkg = "clustrd"),
      varsellcm = list(name = "VarSelLCM", pkg = "VarSelLCM"),
      flexmix = list(name = "FlexMix", pkg = "flexmix")
    )

    valid_methods <- intersect(names(method_meta), selected_methods)
    n_total <- max(1, length(valid_methods))

    # Initialize progress tracker
    init_st <- list()
    for (m in valid_methods) {
      init_st[[m]] <- list(state = "queued", label = method_meta[[m]]$name, detail = "Queued...")
    }
    status_tracker(init_st)

    withProgress(
      message = "Running Fixed-K Horse-Race",
      detail = "Initializing...",
      value = 0,
      {
        idx <- 0
        for (m in valid_methods) {
          idx <- idx + 1
          m_info <- method_meta[[m]]
          m_name <- m_info$name
          m_pkg <- m_info$pkg

          incProgress(
            1 / n_total,
            detail = sprintf("[%d/%d] Running %s (%s)...", idx, n_total, m_name, m_pkg)
          )

          cur_st <- status_tracker()
          cur_st[[m]] <- list(state = "running", label = m_name, detail = "Computing...")
          status_tracker(cur_st)

          if (!has_pkg(m_pkg)) {
            results[[length(results) + 1]] <- data.frame(
              Method = m_name, Package = m_pkg, Time_ms = NA, ARI = NA, Error_Rate = NA,
              Status = "Package Not Installed", stringsAsFactors = FALSE
            )
            cur_st[[m]] <- list(state = "queued", label = m_name, detail = "Not Installed")
            status_tracker(cur_st)
            next
          }

          t_start <- proc.time()
          res_row <- tryCatch({
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
            } else if (m == "cluspcamix") {
              clus_fit <- clustrd::cluspcamix(data = dat$fullData, nclus = k, ndim = 2, nstart = 3)
              memb <- as.integer(clus_fit$cluster)
            } else if (m == "varsellcm") {
              v_fit <- VarSelLCM::VarSelCluster(
                x = dat$fullData, gvals = k, vbleSelec = FALSE, crit.varsel = "BIC", nbcores = 1
              )
              memb <- as.integer(VarSelLCM::fitted(v_fit, type = "partition"))
            } else if (m == "flexmix") {
              cat_mat <- model.matrix(~ . - 1, data = dat$catVars)
              con_mat <- as.matrix(dat$conVars)
              f_fit <- flexmix::flexmix(
                cbind(con_mat, cat_mat) ~ 1,
                k = k,
                model = list(
                  flexmix::FLXMCmvnorm(con_mat ~ 1, diagonal = TRUE),
                  flexmix::FLXMCmvbinary(cat_mat ~ 1)
                ),
                control = list(iter.max = 30, minprior = 0.05, verbose = 0)
              )
              memb <- as.integer(flexmix::clusters(f_fit))
            }

            t_elapsed <- (proc.time() - t_start)[["elapsed"]] * 1000
            ari <- calc_ari(dat$trueID, memb)
            err <- calc_misclass_error(dat$trueID, memb)

            cur_st[[m]] <- list(
              state = "done", label = m_name,
              detail = sprintf("Done (%0.0f ms, ARI = %0.2f)", t_elapsed, ari)
            )
            status_tracker(cur_st)

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
            cur_st[[m]] <- list(state = "error", label = m_name, detail = paste("Error:", e$message))
            status_tracker(cur_st)

            data.frame(
              Method = m_name, Package = m_pkg, Time_ms = NA, ARI = NA, Error_Rate = NA,
              Status = paste("Error:", substr(e$message, 1, 28)), stringsAsFactors = FALSE
            )
          })

          results[[length(results) + 1]] <- res_row
        }
      }
    )

    if (length(results) == 0) {
      return(data.frame(Message = "No techniques selected"))
    }
    do.call(rbind, results)
  }, ignoreNULL = FALSE)

  # ----------------------------------------------------------------------------
  # Cluster Selection Benchmark (K in 2:5)
  # ----------------------------------------------------------------------------
  selection_results <- eventReactive(input$btn_run_select, {
    dat <- sim_data()
    selected_methods <- input$methods
    true_k <- isolate(input$k_clusters)
    results <- list()

    withProgress(
      message = "Running Cluster Selection (K in 2:5)",
      detail = "Testing cluster criteria...",
      value = 0,
      {
        # 1. KAMILA Cluster Selection
        if ("kamila" %in% selected_methods) {
          incProgress(0.16, detail = "Running KAMILA Prediction Strength...")
          if (has_pkg("kamila")) {
            t_start <- proc.time()
            res_k <- tryCatch({
              kamila::kamila(
                dat$conVars,
                dat$catVars,
                numClust = 2:5,
                numInit = 3,
                calcNumClust = "ps",
                numPredStrCvRun = 5,
                predStrThresh = 0.6
              )
            }, error = function(e) NULL)
            t_elapsed <- (proc.time() - t_start)[["elapsed"]] * 1000

            if (!is.null(res_k)) {
              pred_k <- if (is.list(res_k$nClust)) res_k$nClust$bestNClust else res_k$nClust
              ari <- calc_ari(dat$trueID, res_k$finalMemb)
              results[[length(results) + 1]] <- data.frame(
                Method = "KAMILA",
                Package = "kamila",
                True_K = true_k,
                Predicted_K = pred_k,
                Criterion = "Prediction Strength",
                ARI = round(ari, 4),
                Time_ms = round(t_elapsed, 1),
                stringsAsFactors = FALSE
              )
            }
          }
        }

        # 2. Gower + PAM Silhouette Selection
        if ("gower_pam" %in% selected_methods) {
          incProgress(0.16, detail = "Evaluating Gower + PAM Silhouette Widths...")
          if (has_pkg("cluster")) {
            t_start <- proc.time()
            res_pam <- tryCatch({
              g_dist <- cluster::daisy(dat$fullData, metric = "gower")
              sils <- sapply(2:5, function(ki) {
                cluster::pam(g_dist, k = ki, diss = TRUE)$silinfo$avg.width
              })
              best_ki <- (2:5)[which.max(sils)]
              fit_best <- cluster::pam(g_dist, k = best_ki, diss = TRUE)
              list(k = best_ki, memb = as.integer(fit_best$clustering))
            }, error = function(e) NULL)
            t_elapsed <- (proc.time() - t_start)[["elapsed"]] * 1000

            if (!is.null(res_pam)) {
              ari <- calc_ari(dat$trueID, res_pam$memb)
              results[[length(results) + 1]] <- data.frame(
                Method = "Gower + PAM",
                Package = "cluster",
                True_K = true_k,
                Predicted_K = res_pam$k,
                Criterion = "Avg Silhouette Width",
                ARI = round(ari, 4),
                Time_ms = round(t_elapsed, 1),
                stringsAsFactors = FALSE
              )
            }
          }
        }

        # 3. K-Prototypes Validation Index Selection
        if ("kproto" %in% selected_methods) {
          incProgress(0.16, detail = "Running K-Prototypes Validation Indices...")
          if (has_pkg("clustMixType")) {
            t_start <- proc.time()
            res_kp <- tryCatch({
              val <- clustMixType::validation_kproto(
                method = "silhouette",
                data = dat$fullData,
                k = 2:5,
                nstart = 2,
                verbose = FALSE
              )
              best_k <- val$k_opt
              fit_kp <- clustMixType::kproto(dat$fullData, k = best_k, nstart = 2, verbose = FALSE)
              list(k = best_k, memb = as.integer(fit_kp$cluster))
            }, error = function(e) NULL)
            t_elapsed <- (proc.time() - t_start)[["elapsed"]] * 1000

            if (!is.null(res_kp)) {
              ari <- calc_ari(dat$trueID, res_kp$memb)
              results[[length(results) + 1]] <- data.frame(
                Method = "K-Prototypes",
                Package = "clustMixType",
                True_K = true_k,
                Predicted_K = res_kp$k,
                Criterion = "Silhouette Index",
                ARI = round(ari, 4),
                Time_ms = round(t_elapsed, 1),
                stringsAsFactors = FALSE
              )
            }
          }
        }

        # 4. ClusPCAMix Selection
        if ("cluspcamix" %in% selected_methods) {
          incProgress(0.16, detail = "Running ClusPCAMix Selection...")
          if (has_pkg("clustrd")) {
            t_start <- proc.time()
            res_c <- tryCatch({
              fits <- lapply(2:5, function(ki) {
                clustrd::cluspcamix(data = dat$fullData, nclus = ki, ndim = 2, nstart = 2)
              })
              objs <- sapply(fits, function(f) f$criterion)
              best_ki <- (2:5)[which.max(objs)]
              list(k = best_ki, memb = as.integer(fits[[which.max(objs)]]$cluster))
            }, error = function(e) NULL)
            t_elapsed <- (proc.time() - t_start)[["elapsed"]] * 1000

            if (!is.null(res_c)) {
              ari <- calc_ari(dat$trueID, res_c$memb)
              results[[length(results) + 1]] <- data.frame(
                Method = "ClusPCAMix",
                Package = "clustrd",
                True_K = true_k,
                Predicted_K = res_c$k,
                Criterion = "Objective Value",
                ARI = round(ari, 4),
                Time_ms = round(t_elapsed, 1),
                stringsAsFactors = FALSE
              )
            }
          }
        }

        # 5. VarSelLCM Selection
        if ("varsellcm" %in% selected_methods) {
          incProgress(0.16, detail = "Running VarSelLCM Selection (BIC)...")
          if (has_pkg("VarSelLCM")) {
            t_start <- proc.time()
            res_v <- tryCatch({
              VarSelLCM::VarSelCluster(
                x = dat$fullData, gvals = 2:5, vbleSelec = FALSE, crit.varsel = "BIC", nbcores = 1
              )
            }, error = function(e) NULL)
            t_elapsed <- (proc.time() - t_start)[["elapsed"]] * 1000

            if (!is.null(res_v)) {
              memb_v <- as.integer(VarSelLCM::fitted(res_v, type = "partition"))
              pred_k <- tryCatch({
                if (methods::.hasSlot(res_v, "model") && methods::.hasSlot(res_v@model, "g")) {
                  as.integer(res_v@model@g)
                } else {
                  length(unique(memb_v))
                }
              }, error = function(e) length(unique(memb_v)))
              ari <- calc_ari(dat$trueID, memb_v)
              results[[length(results) + 1]] <- data.frame(
                Method = "VarSelLCM",
                Package = "VarSelLCM",
                True_K = true_k,
                Predicted_K = pred_k,
                Criterion = "BIC / MICL",
                ARI = round(ari, 4),
                Time_ms = round(t_elapsed, 1),
                stringsAsFactors = FALSE
              )
            }
          }
        }

        # 6. FlexMix Selection
        if ("flexmix" %in% selected_methods) {
          incProgress(0.16, detail = "Running FlexMix BIC Selection...")
          if (has_pkg("flexmix")) {
            t_start <- proc.time()
            res_f <- tryCatch({
              cat_mat <- model.matrix(~ . - 1, data = dat$catVars)
              con_mat <- as.matrix(dat$conVars)
              m_step <- flexmix::stepFlexmix(
                cbind(con_mat, cat_mat) ~ 1,
                k = 2:5,
                nrep = 2,
                model = list(
                  flexmix::FLXMCmvnorm(con_mat ~ 1, diagonal = TRUE),
                  flexmix::FLXMCmvbinary(cat_mat ~ 1)
                ),
                control = list(iter.max = 20, minprior = 0.05, verbose = 0)
              )
              best_m <- flexmix::getModel(m_step, which = "BIC")
              list(k = best_m@k, memb = as.integer(flexmix::clusters(best_m)))
            }, error = function(e) NULL)
            t_elapsed <- (proc.time() - t_start)[["elapsed"]] * 1000

            if (!is.null(res_f)) {
              ari <- calc_ari(dat$trueID, res_f$memb)
              results[[length(results) + 1]] <- data.frame(
                Method = "FlexMix",
                Package = "flexmix",
                True_K = true_k,
                Predicted_K = res_f$k,
                Criterion = "BIC",
                ARI = round(ari, 4),
                Time_ms = round(t_elapsed, 1),
                stringsAsFactors = FALSE
              )
            }
          }
        }
      }
    )

    if (length(results) == 0) {
      return(data.frame(Message = "Run selection by clicking 'Run Cluster Selection (K in 2:5)'"))
    }
    do.call(rbind, results)
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

    # ClusPCAMix
    if ("cluspcamix" %in% selected_methods && has_pkg("clustrd")) {
      clus_res <- tryCatch({
        as.integer(clustrd::cluspcamix(data = dat$fullData, nclus = k, ndim = 2, nstart = 2)$cluster)
      }, error = function(e) NULL)
      if (!is.null(clus_res)) membs$cluspcamix <- clus_res
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

    # FlexMix
    if ("flexmix" %in% selected_methods && has_pkg("flexmix")) {
      f_res <- tryCatch({
        cat_mat <- model.matrix(~ . - 1, data = dat$catVars)
        con_mat <- as.matrix(dat$conVars)
        m <- flexmix::flexmix(
          cbind(con_mat, cat_mat) ~ 1,
          k = k,
          model = list(
            flexmix::FLXMCmvnorm(con_mat ~ 1, diagonal = TRUE),
            flexmix::FLXMCmvbinary(cat_mat ~ 1)
          ),
          control = list(iter.max = 20, minprior = 0.05, verbose = 0)
        )
        as.integer(flexmix::clusters(m))
      }, error = function(e) NULL)
      if (!is.null(f_res)) membs$flexmix <- f_res
    }

    membs
  })

  # Method code to display name mapping
  method_name_map <- c(
    kamila = "KAMILA",
    gower_pam = "Gower + PAM",
    kproto = "K-Prototypes",
    cluspcamix = "ClusPCAMix",
    varsellcm = "VarSelLCM",
    flexmix = "FlexMix"
  )

  # Render Tables & Plots
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

    par(mfrow = c(1, 2), mar = c(7.5, 4.5, 3, 1))

    # ARI Plot
    barplot(
      valid_res$ARI,
      names.arg = valid_res$Method,
      col = "#2c3e50",
      main = "Adjusted Rand Index (Higher = Better)",
      ylab = "ARI Score",
      ylim = c(0, 1),
      las = 2,
      cex.names = 0.95
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
      cex.names = 0.95
    )
    abline(h = axTicks(2), col = "gray80", lty = 2)
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

    par(mar = c(7.5, 4.5, 3, 1))
    barplot(
      valid_res$Predicted_K,
      names.arg = valid_res$Method,
      col = "#3498db",
      main = "Predicted Number of Clusters (True K indicated by dashed line)",
      ylab = "Predicted K",
      ylim = c(0, 6),
      las = 2,
      cex.names = 0.95
    )
    abline(h = isolate(input$k_clusters), col = "red", lty = 2, lwd = 2)
    legend("topright", legend = paste("True K =", isolate(input$k_clusters)), col = "red", lty = 2, lwd = 2)
  })

  # ----------------------------------------------------------------------------
  # LDA Cluster Projection Plot (True Points LDA with Multi-Method Assignments)
  # ----------------------------------------------------------------------------
  output$lda_cluster_plot <- renderPlot({
    dat <- sim_data()
    all_membs <- cluster_assignments()
    k_true <- isolate(input$k_clusters)

    # Construct design matrix for true LDA projection
    cat_mat <- model.matrix(~ ., data = dat$catVars)[, -1, drop = FALSE]
    comb_mat <- cbind(scale(as.matrix(dat$conVars)), scale(cat_mat))
    zv <- apply(comb_mat, 2, function(x) var(x, na.rm = TRUE) == 0 || is.na(var(x)))
    comb_mat <- comb_mat[, !zv, drop = FALSE]

    df_lda <- as.data.frame(comb_mat)
    df_lda$class <- as.factor(dat$trueID)

    # Compute LDA using true ground truth classes
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
      coords_y <- dat$conVars[, 1]
      xlab_txt <- "Linear Discriminant 1 (LD1)"
      ylab_txt <- "Continuous Feature 1 (Con_1)"
    }

    palette <- c("#e74c3c", "#3498db", "#2ecc71", "#f39c12", "#9b59b6", "#1abc9c", "#e67e22")

    method_labels <- c(
      true = "True Ground Truth",
      kamila = "KAMILA (kamila)",
      gower_pam = "Gower + PAM (cluster)",
      kproto = "K-Prototypes (clustMixType)",
      cluspcamix = "ClusPCAMix (clustrd)",
      varsellcm = "VarSelLCM",
      flexmix = "FlexMix"
    )

    if (input$proj_view == "single") {
      # Single Focused Plot
      target_m <- input$proj_single_method
      cur_memb <- if (target_m %in% names(all_membs)) all_membs[[target_m]] else dat$trueID
      title_str <- method_labels[target_m]
      if (is.na(title_str)) title_str <- target_m

      par(mar = c(5, 5, 4, 2))
      plot(
        coords_x, coords_y,
        col = palette[cur_memb],
        pch = 19,
        cex = 1.0,
        main = paste0("LDA Discriminant Projection: ", title_str, " (N = ", length(coords_x), ")"),
        xlab = xlab_txt,
        ylab = ylab_txt
      )
      grid()
      legend(
        "topright",
        legend = paste("Cluster", seq_len(k_true)),
        col = palette[seq_len(k_true)],
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
        cur_memb <- all_membs[[m_key]]
        m_title <- method_labels[m_key]
        if (is.na(m_title)) m_title <- m_key

        plot(
          coords_x, coords_y,
          col = palette[cur_memb],
          pch = if (m_key == "true") 19 else 17,
          cex = 0.8,
          main = m_title,
          xlab = xlab_txt,
          ylab = ylab_txt
        )
        grid()
      }
    }
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
    pkgs <- c("kamila", "cluster", "clustMixType", "clustrd", "VarSelLCM", "flexmix", "mixtools", "mclust", "shiny")
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
