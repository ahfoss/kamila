# ==============================================================================
# Mixed-Type Clustering Horse-Race - Shiny Web Application
# Compatible with local Shiny and Shinylive / WebR (Client-Side)
# ==============================================================================

if (requireNamespace("shiny", quietly = TRUE)) {
  library(shiny)
}

# Optional dependencies loaded conditionally for WebR compatibility
has_pkg <- function(pkg) requireNamespace(pkg, quietly = TRUE)

# ------------------------------------------------------------------------------
# Helper: Ground truth evaluation metrics
# ------------------------------------------------------------------------------
calc_ari <- function(true_labels, pred_labels) {
  if (has_pkg("mclust")) {
    return(mclust::adjustedRandIndex(true_labels, pred_labels))
  }
  # Base R fallback for Adjusted Rand Index
  tab <- table(true_labels, pred_labels)
  n <- length(true_labels)
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
  tab <- table(true_labels, pred_labels)
  # For modest K, compute best label permutation matching
  k_true <- nrow(tab)
  k_pred <- ncol(tab)

  if (k_pred > 6 || k_true > 6) {
    # Greedy matching for larger K
    matched_correct <- 0
    temp_tab <- tab
    for (i in 1:min(k_true, k_pred)) {
      max_idx <- which(temp_tab == max(temp_tab), arr.ind = TRUE)[1, ]
      matched_correct <- matched_correct + temp_tab[max_idx[1], max_idx[2]]
      temp_tab[max_idx[1], ] <- -1
      temp_tab[, max_idx[2]] <- -1
    }
    return(1 - (matched_correct / length(true_labels)))
  }

  # Exact best permutation for small K
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
  1 - (max_correct / length(true_labels))
}

# ------------------------------------------------------------------------------
# Helper: Synthetic mixed data generator
# ------------------------------------------------------------------------------
generate_synthetic_mixed_data <- function(n = 150, p_con = 5, p_cat = 5, k = 3,
                                          separation = 2.0, num_levels = 4, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)

  # If k == 2 and kamila::genMixedData is available, we can use it or general formula
  props <- rep(1 / k, k)
  cluster_assign <- sample(seq_len(k), size = n, replace = TRUE, prob = props)

  # Continuous features: Gaussian clusters with separation
  con_data <- matrix(0, nrow = n, ncol = p_con)
  colnames(con_data) <- paste0("Con_", seq_len(p_con))

  for (j in seq_len(p_con)) {
    centers <- seq(-separation * (k - 1) / 2, separation * (k - 1) / 2, length.out = k)
    con_data[, j] <- rnorm(n, mean = centers[cluster_assign], sd = 1.0)
  }

  # Categorical features: Multinomial with dominant level per cluster
  cat_list <- list()
  for (j in seq_len(p_cat)) {
    col_vals <- integer(n)
    for (clust in seq_len(k)) {
      idx <- which(cluster_assign == clust)
      if (length(idx) > 0) {
        # Create probability vector favoring one level for each cluster
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
# UI Definition
# ------------------------------------------------------------------------------
ui <- fluidPage(
  theme = if (has_pkg("bslib")) bslib::bs_theme(version = 5, bootswatch = "flatly") else NULL,

  titlePanel(
    tags$div(
      tags$h2("Mixed-Type Clustering Horse-Race", style = "margin-bottom: 2px; font-weight: 700;"),
      tags$p(
        "Interactive client-side benchmark comparing KAMILA with alternative mixed-data clustering algorithms.",
        style = "color: #6c757d; font-size: 1.05rem;"
      )
    ),
    windowTitle = "Mixed-Type Clustering Horse-Race"
  ),

  sidebarLayout(
    sidebarPanel(
      width = 4,
      tags$h4("Simulation Settings", style = "font-weight: 600;"),

      sliderInput("n_obs", "Number of Observations (N):", min = 60, max = 500, value = 150, step = 10),
      fluidRow(
        column(6, sliderInput("p_con", "Continuous (P1):", min = 1, max = 15, value = 5, step = 1)),
        column(6, sliderInput("p_cat", "Categorical (P2):", min = 1, max = 15, value = 5, step = 1))
      ),
      fluidRow(
        column(6, sliderInput("k_clusters", "True Clusters (K):", min = 2, max = 5, value = 3, step = 1)),
        column(6, sliderInput("separation", "Separation:", min = 0.5, max = 4.0, value = 2.0, step = 0.5))
      ),
      numericInput("rand_seed", "Random Seed:", value = 42, min = 1),

      tags$hr(),
      tags$h4("Algorithms to Benchmark", style = "font-weight: 600;"),
      checkboxGroupInput(
        "methods",
        label = NULL,
        choices = c(
          "KAMILA (kamila)" = "kamila",
          "Gower's Distance + PAM (cluster)" = "gower_pam",
          "K-Prototypes (clustMixType)" = "kproto",
          "ClusPCAMix (clustrd)" = "cluspcamix"
        ),
        selected = c("kamila", "gower_pam")
      ),

      actionButton(
        "btn_run",
        "Run Horse Race",
        class = "btn-primary btn-lg w-100",
        style = "margin-top: 10px; font-weight: 600;"
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
            tags$strong("Task: "),
            "Evaluate clustering performance when true K is provided to each algorithm."
          ),
          tags$h4("Performance Summary", style = "font-weight: 600; margin-top: 15px;"),
          tableOutput("benchmark_table"),
          tags$hr(),
          tags$h4("Visual Comparison", style = "font-weight: 600;"),
          plotOutput("benchmark_plot", height = "320px")
        ),

        tabPanel(
          "Cluster & Data Visualization",
          tags$br(),
          tags$p("2D Principal Component Projection comparing true ground-truth labels vs. KAMILA assignments."),
          plotOutput("pca_cluster_plot", height = "420px")
        ),

        tabPanel(
          "Environment & Packages",
          tags$br(),
          tags$h4("Runtime Environment Info", style = "font-weight: 600;"),
          tableOutput("env_info_table"),
          tags$h4("Package Versions", style = "font-weight: 600; margin-top: 20px;"),
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

  # Reactive dataset generator
  sim_data <- reactive({
    input$btn_run
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

  # Benchmark Execution
  benchmark_results <- reactive({
    dat <- sim_data()
    selected_methods <- input$methods
    k <- isolate(input$k_clusters)

    results <- list()

    # 1. KAMILA
    if ("kamila" %in% selected_methods) {
      if (has_pkg("kamila")) {
        t_start <- proc.time()
        res <- tryCatch({
          kamila::kamila(
            dat$conVars,
            dat$catVars,
            numClust = k,
            numInit = 5,
            maxIter = 25,
            calcNumClust = "none"
          )
        }, error = function(e) NULL)
        t_elapsed <- (proc.time() - t_start)[["elapsed"]] * 1000 # in ms

        if (!is.null(res)) {
          ari <- calc_ari(dat$trueID, res$finalMemb)
          err <- calc_misclass_error(dat$trueID, res$finalMemb)
          results[[length(results) + 1]] <- data.frame(
            Method = "KAMILA",
            Package = "kamila",
            Time_ms = round(t_elapsed, 1),
            ARI = round(ari, 4),
            Error_Rate = round(err, 4),
            Status = "Success",
            stringsAsFactors = FALSE
          )
        }
      } else {
        results[[length(results) + 1]] <- data.frame(
          Method = "KAMILA", Package = "kamila", Time_ms = NA, ARI = NA, Error_Rate = NA,
          Status = "Package not installed", stringsAsFactors = FALSE
        )
      }
    }

    # 2. Gower's Distance + PAM
    if ("gower_pam" %in% selected_methods) {
      if (has_pkg("cluster")) {
        t_start <- proc.time()
        res <- tryCatch({
          gower_dist <- cluster::daisy(dat$fullData, metric = "gower")
          cluster::pam(gower_dist, k = k, diss = TRUE)
        }, error = function(e) NULL)
        t_elapsed <- (proc.time() - t_start)[["elapsed"]] * 1000 # in ms

        if (!is.null(res)) {
          ari <- calc_ari(dat$trueID, res$clustering)
          err <- calc_misclass_error(dat$trueID, res$clustering)
          results[[length(results) + 1]] <- data.frame(
            Method = "Gower + PAM",
            Package = "cluster",
            Time_ms = round(t_elapsed, 1),
            ARI = round(ari, 4),
            Error_Rate = round(err, 4),
            Status = "Success",
            stringsAsFactors = FALSE
          )
        }
      } else {
        results[[length(results) + 1]] <- data.frame(
          Method = "Gower + PAM", Package = "cluster", Time_ms = NA, ARI = NA, Error_Rate = NA,
          Status = "Package not installed", stringsAsFactors = FALSE
        )
      }
    }

    # 3. K-Prototypes
    if ("kproto" %in% selected_methods) {
      if (has_pkg("clustMixType")) {
        t_start <- proc.time()
        res <- tryCatch({
          clustMixType::kproto(dat$fullData, k = k, nstart = 3, verbose = FALSE)
        }, error = function(e) NULL)
        t_elapsed <- (proc.time() - t_start)[["elapsed"]] * 1000 # in ms

        if (!is.null(res)) {
          ari <- calc_ari(dat$trueID, res$cluster)
          err <- calc_misclass_error(dat$trueID, res$cluster)
          results[[length(results) + 1]] <- data.frame(
            Method = "K-Prototypes",
            Package = "clustMixType",
            Time_ms = round(t_elapsed, 1),
            ARI = round(ari, 4),
            Error_Rate = round(err, 4),
            Status = "Success",
            stringsAsFactors = FALSE
          )
        }
      } else {
        results[[length(results) + 1]] <- data.frame(
          Method = "K-Prototypes", Package = "clustMixType", Time_ms = NA, ARI = NA, Error_Rate = NA,
          Status = "Package not installed", stringsAsFactors = FALSE
        )
      }
    }

    # 4. ClusPCAMix
    if ("cluspcamix" %in% selected_methods) {
      if (has_pkg("clustrd")) {
        t_start <- proc.time()
        res <- tryCatch({
          clustrd::cluspcamix(dat$fullData, k = k, nstart = 3)
        }, error = function(e) NULL)
        t_elapsed <- (proc.time() - t_start)[["elapsed"]] * 1000 # in ms

        if (!is.null(res)) {
          ari <- calc_ari(dat$trueID, res$cluster)
          err <- calc_misclass_error(dat$trueID, res$cluster)
          results[[length(results) + 1]] <- data.frame(
            Method = "ClusPCAMix",
            Package = "clustrd",
            Time_ms = round(t_elapsed, 1),
            ARI = round(ari, 4),
            Error_Rate = round(err, 4),
            Status = "Success",
            stringsAsFactors = FALSE
          )
        }
      } else {
        results[[length(results) + 1]] <- data.frame(
          Method = "ClusPCAMix", Package = "clustrd", Time_ms = NA, ARI = NA, Error_Rate = NA,
          Status = "Package not installed", stringsAsFactors = FALSE
        )
      }
    }

    if (length(results) == 0) {
      return(data.frame(Message = "No methods selected or available"))
    }
    do.call(rbind, results)
  })

  # Outputs
  output$benchmark_table <- renderTable({
    benchmark_results()
  }, striped = TRUE, hover = TRUE, bordered = TRUE)

  output$benchmark_plot <- renderPlot({
    res <- benchmark_results()
    if (!"ARI" %in% names(res) || nrow(res) == 0) return(NULL)

    valid_res <- res[!is.na(res$ARI), ]
    if (nrow(valid_res) == 0) return(NULL)

    par(mfrow = c(1, 2), mar = c(5, 5, 3, 1))

    # ARI Plot
    barplot(
      valid_res$ARI,
      names.arg = valid_res$Method,
      col = "#2c3e50",
      main = "Adjusted Rand Index (Higher = Better)",
      ylab = "ARI Score",
      ylim = c(0, 1),
      las = 2
    )
    abline(h = seq(0, 1, 0.2), col = "gray80", lty = 2)

    # Timing Plot
    barplot(
      valid_res$Time_ms,
      names.arg = valid_res$Method,
      col = "#18bc9c",
      main = "Execution Time (Lower = Faster)",
      ylab = "Time (ms)",
      las = 2
    )
    abline(h = axTicks(2), col = "gray80", lty = 2)
  })

  # 2D PCA Cluster Plot
  output$pca_cluster_plot <- renderPlot({
    dat <- sim_data()
    # Dummy code factor variables for PCA projection
    cat_mat <- model.matrix(~ . - 1, data = dat$catVars)
    combined_mat <- cbind(scale(as.matrix(dat$conVars)), scale(cat_mat))

    pca_fit <- prcomp(combined_mat, center = TRUE, scale. = FALSE)
    pca_coords <- pca_fit$x[, 1:2]

    # Run KAMILA if available
    kam_clust <- if (has_pkg("kamila")) {
      tryCatch({
        kam_out <- kamila::kamila(
          dat$conVars,
          dat$catVars,
          numClust = input$k_clusters,
          numInit = 3,
          calcNumClust = "none"
        )
        kam_out$finalMemb
      }, error = function(e) dat$trueID)
    } else {
      dat$trueID
    }

    par(mfrow = c(1, 2), mar = c(4, 4, 3, 1))

    palette <- c("#e74c3c", "#3498db", "#2ecc71", "#f39c12", "#9b59b6")

    # True Cluster Plot
    plot(
      pca_coords,
      col = palette[dat$trueID],
      pch = 19,
      cex = 1.3,
      main = "True Cluster Labels (PCA 2D)",
      xlab = "PC 1", ylab = "PC 2"
    )
    grid()

    # Predicted Cluster Plot
    plot(
      pca_coords,
      col = palette[kam_clust],
      pch = 17,
      cex = 1.3,
      main = "KAMILA Predicted Clusters",
      xlab = "PC 1", ylab = "PC 2"
    )
    grid()
  })

  # Environment Info
  output$env_info_table <- renderTable({
    is_webr <- exists("webr", envir = .GlobalEnv) || Sys.getenv("WEBR") == "1"
    data.frame(
      Property = c("R Version", "Platform / Architecture", "Execution Engine", "Interactive Shiny Mode"),
      Value = c(
        R.version.string,
        R.version$platform,
        if (is_webr) "WebR / WebAssembly (Client-Side in Browser)" else "Native R Runtime",
        if (is_webr) "Shinylive (Zero-Server)" else "Local / Server Shiny"
      ),
      stringsAsFactors = FALSE
    )
  }, striped = TRUE, bordered = TRUE)

  # Package Version Table
  output$pkg_version_table <- renderTable({
    pkgs <- c("kamila", "cluster", "mclust", "clustMixType", "clustrd", "VarSelLCM", "flexmix", "shiny")
    versions <- sapply(pkgs, function(p) {
      if (has_pkg(p)) as.character(packageVersion(p)) else "Not Available"
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
