###############################################################
# POLYA_REVISION_PIPELINE.R  (SELF-CONTAINED + COLORFUL OUTPUTS)
# -------------------------------------------------------------
# Outputs auto-saved to:
#   outputs_polya_revision/run_YYYYMMDD_HHMMSS/...
# and copied to stable:
#   outputs_polya_revision/latest/...
#
# Reviewers' requirements satisfied by outputs:
#  - Table with Bias + SE (LaTeX tabular wrapped in revblue)
#  - MovieLens: timestamp-ordered sequential evaluation (LPS + diagnostics)
#  - Dune: current-practice-compatible Dirichlet posterior mean (closed form)
###############################################################

############################
# 0) USER CONFIG
############################
CONFIG <- list(
  seed = 20260223,
  out_base = "outputs_polya_revision",
  
  polya = list(
    a = 10L, b = 10L, c = 1L,
    n_steps = 1000L,
    R = 1000L
  ),
  
  movielens = list(
    local_dir = NULL,        # if you have ml-100k folder, set to that folder path
    train_frac = 0.80,
    alpha0 = 1.0,
    make_plots = TRUE,
    zip_url = "https://files.grouplens.org/datasets/movielens/ml-100k.zip"
  ),
  
  dune = list(
    alpha0 = 1.0,
    make_plots = TRUE
  )
)

############################
# 1) HELPERS
############################
ts_string <- function() format(Sys.time(), "%Y%m%d_%H%M%S")

ensure_dir <- function(path) {
  if (!dir.exists(path)) dir.create(path, recursive = TRUE, showWarnings = FALSE)
  invisible(path)
}

fmt_num <- function(x, digits = 4) formatC(x, format = "f", digits = digits)

write_csv <- function(df, path) {
  utils::write.csv(df, path, row.names = FALSE)
  invisible(path)
}

write_lines <- function(lines, path) {
  writeLines(lines, con = path)
  invisible(path)
}

safe_install <- function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    message("Installing package: ", pkg)
    install.packages(pkg, repos = "https://cloud.r-project.org")
  }
}

# Nice color palettes (base R)
pal_set <- function(n, name = "Set2") grDevices::hcl.colors(n, palette = name)

# Robust one-way counts
count_df <- function(x, name = "value") {
  tab <- table(x, useNA = "ifany")
  df <- data.frame(value = names(tab), count = as.integer(tab), stringsAsFactors = FALSE)
  names(df)[1] <- name
  df
}

# LaTeX tabular wrapped in revblue (for your blue revision workflow)
write_latex_revblue_tabular <- function(df, file, digits = 4) {
  stopifnot(all(c("Parameter", "True", "Mean", "Bias", "SE") %in% names(df)))
  
  lines <- c(
    "\\begin{revblue}",
    "\\begin{tabular}{lrrrr}",
    "\\toprule",
    "Parameter & True & Mean & Bias & SE \\\\",
    "\\midrule"
  )
  
  for (i in seq_len(nrow(df))) {
    row <- df[i, ]
    lines <- c(lines, paste0(
      row$Parameter, " & ",
      fmt_num(row$True, digits), " & ",
      fmt_num(row$Mean, digits), " & ",
      fmt_num(row$Bias, digits), " & ",
      fmt_num(row$SE, digits), " \\\\"
    ))
  }
  
  lines <- c(lines,
             "\\bottomrule",
             "\\end{tabular}",
             "\\end{revblue}"
  )
  write_lines(lines, file)
}

# Copy everything from run_dir -> latest (stable LaTeX paths)
sync_latest <- function(run_root, latest_root) {
  if (dir.exists(latest_root)) unlink(latest_root, recursive = TRUE, force = TRUE)
  dir.create(latest_root, recursive = TRUE, showWarnings = FALSE)
  # Copy full folder
  file.copy(from = run_root, to = latest_root, recursive = TRUE)
  invisible(TRUE)
}

############################
# 2) (A) CLASSICAL PÓLYA URN (Monte Carlo)
############################
polya_one_run <- function(a, b, c, n_steps) {
  # IMPORTANT: Do NOT use NA/NB variable names (NA is reserved in R)
  nA <- as.integer(a)
  nB <- as.integer(b)
  
  for (t in seq_len(n_steps)) {
    pA <- nA / (nA + nB)
    if (runif(1) < pA) nA <- nA + c else nB <- nB + c
  }
  nA / (nA + nB)
}

run_polya_mc <- function(a, b, c, n_steps, R, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  
  est <- numeric(R)
  for (r in seq_len(R)) est[r] <- polya_one_run(a, b, c, n_steps)
  
  # Truthful: E[p_A(n)] = a/(a+b) for all n (martingale)
  true_mean <- a / (a + b)
  
  mean_hat <- mean(est)
  sd_hat   <- sd(est)
  bias     <- mean_hat - true_mean
  se       <- sd_hat / sqrt(R)
  
  summary_df <- data.frame(
    Parameter = sprintf("p_A(%d) (MC mean)", n_steps),
    True = true_mean,
    Mean = mean_hat,
    Bias = bias,
    SE = se,
    check.names = FALSE
  )
  
  list(estimates = est, summary = summary_df)
}

save_polya_outputs <- function(mc, cfg, dirs) {
  # data
  write_csv(data.frame(rep = seq_along(mc$estimates), pA_final = mc$estimates),
            file.path(dirs$data, "polya_mc_estimates.csv"))
  
  # tables
  write_csv(mc$summary, file.path(dirs$tables, "polya_mc_summary_bias_se.csv"))
  
  # latex
  write_latex_revblue_tabular(mc$summary, file.path(dirs$latex, "table_mc_bias_se_revblue.tex"), digits = 4)
  
  # colorful histogram + beta overlay (limit Beta(a/c,b/c))
  png(file.path(dirs$figs, "polya_final_pA_hist_beta_overlay.png"),
      width = 1400, height = 850, res = 160)
  op <- par(mar = c(5,5,4,2) + 0.1)
  cols <- pal_set(3, "Set2")
  hist(mc$estimates, breaks = 40, freq = FALSE,
       col = cols[1], border = "white",
       main = sprintf("Classical Pólya urn: final p_A after n=%d steps (R=%d)", cfg$n_steps, cfg$R),
       xlab = expression(p[A](n)), ylab = "Density")
  # overlay Beta density (limit)
  x <- seq(0, 1, length.out = 500)
  shape1 <- cfg$a / cfg$c
  shape2 <- cfg$b / cfg$c
  lines(x, dbeta(x, shape1 = shape1, shape2 = shape2), lwd = 3, col = cols[2])
  abline(v = cfg$a/(cfg$a+cfg$b), lty = 2, lwd = 2, col = cols[3])
  legend("topright",
         legend = c("MC histogram (density)", sprintf("Beta(%g,%g) limit", shape1, shape2), "True mean a/(a+b)"),
         lwd = c(10,3,2), lty = c(1,1,2), col = c(cols[1], cols[2], cols[3]), bty = "n")
  par(op)
  dev.off()
  
  # log
  txt <- c(
    "Classical Pólya Urn Monte Carlo Summary",
    "--------------------------------------",
    sprintf("a=%d, b=%d, c=%d, n_steps=%d, R=%d", cfg$a, cfg$b, cfg$c, cfg$n_steps, cfg$R),
    sprintf("True E[p_A(n)] = a/(a+b) = %.6f", cfg$a/(cfg$a+cfg$b)),
    sprintf("MC Mean = %.6f", mc$summary$Mean),
    sprintf("Bias    = %.6f", mc$summary$Bias),
    sprintf("SE      = %.6f", mc$summary$SE),
    "",
    "LaTeX table: latex/table_mc_bias_se_revblue.tex",
    "Figure: figs/polya_final_pA_hist_beta_overlay.png"
  )
  write_lines(txt, file.path(dirs$logs, "polya_mc_summary.txt"))
}

############################
# 3) (B) MOVIELENS (timestamps): auto-download + sequential evaluation
############################
download_movielens_100k <- function(zip_url, dest_dir) {
  ensure_dir(dest_dir)
  zip_path <- file.path(dest_dir, "ml-100k.zip")
  ml_dir   <- file.path(dest_dir, "ml-100k")
  udata    <- file.path(ml_dir, "u.data")
  
  if (file.exists(udata)) return(udata)
  
  message("Downloading MovieLens 100k zip to: ", zip_path)
  utils::download.file(zip_url, destfile = zip_path, mode = "wb", quiet = TRUE)
  
  message("Unzipping MovieLens into: ", dest_dir)
  utils::unzip(zip_path, exdir = dest_dir)
  
  if (!file.exists(udata)) stop("After unzip, u.data not found at: ", udata)
  udata
}

run_movielens_seq_polya <- function(u_data_path, train_frac = 0.8, alpha0 = 1.0, make_plots = TRUE, dirs = NULL) {
  dat <- read.table(u_data_path, header = FALSE, stringsAsFactors = FALSE)
  colnames(dat) <- c("user", "item", "rating", "timestamp")
  
  # chronological order
  dat <- dat[order(dat$timestamp, dat$user, dat$item), ]
  
  items <- as.integer(dat$item)
  n <- nrow(dat)
  M <- max(items)
  
  n_train <- floor(train_frac * n)
  train_items <- items[1:n_train]
  test_items  <- items[(n_train + 1):n]
  Ttest <- length(test_items)
  
  train_counts <- tabulate(train_items, nbins = M)
  
  # sequential state + smoothing
  C <- train_counts + alpha0
  sumC <- sum(C)
  
  # baseline fixed multinomial from training
  p0 <- (train_counts + alpha0) / sum(train_counts + alpha0)
  
  logp_polya <- numeric(Ttest)
  logp_base  <- numeric(Ttest)
  
  E_polya <- numeric(M)
  E_base  <- p0 * Ttest
  
  for (t in seq_len(Ttest)) {
    i <- test_items[t]
    p_vec <- C / sumC
    logp_polya[t] <- log(p_vec[i])
    logp_base[t]  <- log(p0[i])
    E_polya <- E_polya + p_vec
    
    # reinforcement
    C[i] <- C[i] + 1
    sumC <- sumC + 1
  }
  
  LPS_polya <- mean(logp_polya)
  LPS_base  <- mean(logp_base)
  
  O <- tabulate(test_items, nbins = M)
  
  chisq_polya_all <- sum((O - E_polya)^2 / E_polya)
  chisq_base_all  <- sum((O - E_base )^2 / E_base )
  
  keep_polya <- which(E_polya >= 5)
  keep_base  <- which(E_base  >= 5)
  
  chisq_polya_ge5 <- sum((O[keep_polya] - E_polya[keep_polya])^2 / E_polya[keep_polya])
  chisq_base_ge5  <- sum((O[keep_base ] - E_base [keep_base ])^2 / E_base [keep_base ])
  
  df_polya_ge5 <- length(keep_polya) - 1
  df_base_ge5  <- length(keep_base ) - 1
  
  diag_df <- data.frame(
    MovieID = seq_len(M),
    Observed = O,
    Expected_Polya = E_polya,
    Expected_Base  = E_base
  )
  
  if (!is.null(dirs)) {
    write_csv(diag_df, file.path(dirs$tables, "movielens_observed_expected_counts.csv"))
    
    # LaTeX summary table
    latex_lines <- c(
      "\\begin{tabular}{lr}",
      "\\toprule",
      "Quantity & Value \\\\",
      "\\midrule",
      paste0("Total events $n$ & ", n, " \\\\"),
      paste0("Train events & ", n_train, " \\\\"),
      paste0("Test events $T$ & ", Ttest, " \\\\"),
      paste0("Movies $M$ & ", M, " \\\\"),
      paste0("Smoothing $\\alpha_0$ & ", fmt_num(alpha0, 3), " \\\\"),
      paste0("LPS (Sequential P\\'olya) & ", fmt_num(LPS_polya, 6), " \\\\"),
      paste0("LPS (Baseline) & ", fmt_num(LPS_base, 6), " \\\\"),
      paste0("$\\chi^2$ P\\'olya (all bins) & ", fmt_num(chisq_polya_all, 4), " \\\\"),
      paste0("$\\chi^2$ Baseline (all bins) & ", fmt_num(chisq_base_all, 4), " \\\\"),
      paste0("$\\chi^2$ P\\'olya ($E_i\\ge 5$) & ", fmt_num(chisq_polya_ge5, 4), " \\\\"),
      paste0("df (P\\'olya; $E_i\\ge 5$) & ", df_polya_ge5, " \\\\"),
      paste0("$\\chi^2$ Baseline ($E_i\\ge 5$) & ", fmt_num(chisq_base_ge5, 4), " \\\\"),
      paste0("df (Baseline; $E_i\\ge 5$) & ", df_base_ge5, " \\\\"),
      "\\bottomrule",
      "\\end{tabular}"
    )
    write_lines(latex_lines, file.path(dirs$latex, "table_movielens_sequential_eval.tex"))
    
    # colorful plots
    if (make_plots) {
      cols <- pal_set(3, "Dark 3")
      
      png(file.path(dirs$figs, "movielens_running_mean_logpred.png"),
          width = 1400, height = 850, res = 160)
      op <- par(mar = c(5,5,4,2) + 0.1)
      rm_polya <- cumsum(logp_polya) / seq_along(logp_polya)
      rm_base  <- cumsum(logp_base)  / seq_along(logp_base)
      plot(rm_polya, type = "l", lwd = 3, col = cols[1],
           xlab = "Test step (chronological)",
           ylab = "Running mean log predictive probability",
           main = "MovieLens: sequential one-step-ahead prediction (timestamp-ordered)")
      lines(rm_base, lwd = 3, col = cols[2], lty = 2)
      legend("bottomright",
             legend = c("Sequential Pólya", "Baseline (fixed train multinomial)"),
             col = c(cols[1], cols[2]), lty = c(1,2), lwd = 3, bty = "n")
      par(op)
      dev.off()
      
      # Top-30 movies observed vs expected (barplot)
      topK <- 30
      ord <- order(O, decreasing = TRUE)
      top_ids <- ord[1:topK]
      mat <- rbind(O[top_ids], E_polya[top_ids], E_base[top_ids])
      rownames(mat) <- c("Observed", "Expected (Pólya)", "Expected (Baseline)")
      
      bar_cols <- pal_set(nrow(mat), "Set2")
      
      png(file.path(dirs$figs, "movielens_top_movies_observed_expected.png"),
          width = 1600, height = 900, res = 160)
      op <- par(mar = c(10,5,4,2) + 0.1)
      bp <- barplot(mat, beside = TRUE, col = bar_cols, border = "white",
                    names.arg = paste0("M", top_ids),
                    las = 2, cex.names = 0.65,
                    main = sprintf("Held-out block: observed vs expected counts (Top %d movies)", topK),
                    ylab = "Count")
      legend("topright", legend = rownames(mat), fill = bar_cols, bty = "n")
      par(op)
      dev.off()
      
      # Scatter observed vs expected (log scale)
      png(file.path(dirs$figs, "movielens_observed_vs_expected_scatter.png"),
          width = 1400, height = 850, res = 160)
      op <- par(mar = c(5,5,4,2) + 0.1)
      x <- log10(1 + diag_df$Expected_Polya)
      y <- log10(1 + diag_df$Observed)
      plot(x, y, pch = 16, cex = 0.6, col = grDevices::adjustcolor(cols[1], 0.35),
           xlab = expression(log[10](1 + E[i]~"(Pólya)")),
           ylab = expression(log[10](1 + O[i])),
           main = "MovieLens: Observed vs Expected counts (Pólya), log-scale")
      abline(0, 1, lty = 2, lwd = 2, col = cols[3])
      par(op)
      dev.off()
    }
    
    # log
    write_lines(c(
      "MovieLens sequential (timestamp) evaluation",
      "------------------------------------------",
      sprintf("u.data: %s", u_data_path),
      sprintf("n_total=%d, n_train=%d, n_test=%d, M=%d", n, n_train, Ttest, M),
      sprintf("alpha0=%.3f, train_frac=%.2f", alpha0, train_frac),
      sprintf("LPS_polya=%.6f", LPS_polya),
      sprintf("LPS_base =%.6f", LPS_base),
      sprintf("chisq_polya_all=%.4f", chisq_polya_all),
      sprintf("chisq_base_all =%.4f", chisq_base_all),
      sprintf("chisq_polya_ge5=%.4f (df=%d)", chisq_polya_ge5, df_polya_ge5),
      sprintf("chisq_base_ge5 =%.4f (df=%d)", chisq_base_ge5, df_base_ge5)
    ), file.path(dirs$logs, "movielens_sequential_eval_summary.txt"))
  }
  
  invisible(list(
    n_total = n, n_train = n_train, n_test = Ttest, n_movies = M,
    alpha0 = alpha0, train_frac = train_frac,
    LPS_polya = LPS_polya, LPS_base = LPS_base,
    chisq_polya_all = chisq_polya_all,
    chisq_base_all  = chisq_base_all,
    chisq_polya_ge5 = chisq_polya_ge5,
    chisq_base_ge5  = chisq_base_ge5,
    df_polya_ge5 = df_polya_ge5,
    df_base_ge5  = df_base_ge5
  ))
}

############################
# 4) (C) MASS::crabs outputs (correct variables: sp, sex)
############################
run_crabs_outputs <- function(dirs) {
  safe_install("MASS")
  data("crabs", package = "MASS")
  cr <- MASS::crabs
  
  write_csv(cr, file.path(dirs$data, "crabs_raw.csv"))
  write_csv(count_df(cr$sp,  "sp"),  file.path(dirs$tables, "crabs_sp_counts.csv"))
  write_csv(count_df(cr$sex, "sex"), file.path(dirs$tables, "crabs_sex_counts.csv"))
  
  # simple Pólya demo on species (illustrative)
  labels <- levels(as.factor(cr$sp))
  init_counts <- tabulate(as.integer(as.factor(cr$sp)), nbins = length(labels))
  
  alpha0 <- 1
  state <- init_counts + alpha0
  T <- nrow(cr)
  draws <- integer(T)
  
  for (t in seq_len(T)) {
    p <- state / sum(state)
    draws[t] <- sample.int(length(state), size = 1, prob = p)
    state[draws[t]] <- state[draws[t]] + 1
  }
  sim_counts <- tabulate(draws, nbins = length(state))
  
  demo <- data.frame(sp = labels, observed = init_counts, polya_demo_simulated = sim_counts)
  write_csv(demo, file.path(dirs$tables, "crabs_sp_polya_demo.csv"))
  
  cols <- pal_set(2, "Set2")
  png(file.path(dirs$figs, "crabs_sp_observed_vs_polya_demo.png"),
      width = 1200, height = 750, res = 160)
  op <- par(mar = c(5,5,4,2) + 0.1)
  mat <- rbind(init_counts, sim_counts)
  colnames(mat) <- labels
  rownames(mat) <- c("Observed", "Pólya demo")
  barplot(mat, beside = TRUE, col = cols, border = "white",
          main = "MASS::crabs: species counts (Observed vs Pólya demo)",
          ylab = "Count")
  legend("topright", legend = rownames(mat), fill = cols, bty = "n")
  par(op)
  dev.off()
}

############################
# 5) (D) vegan::dune Dirichlet posterior mean (closed form)
############################
run_dune_dirichlet_posterior <- function(alpha0 = 1.0, make_plots = TRUE, dirs = NULL) {
  safe_install("vegan")
  data("dune", package = "vegan")
  dune_mat <- dune
  
  counts <- colSums(dune_mat)
  species <- names(counts)
  N <- sum(counts)
  
  theta_true <- counts / N
  alpha_post <- counts + alpha0
  theta_post_mean <- alpha_post / sum(alpha_post)
  
  out <- data.frame(
    Species = species,
    True_Proportion = as.numeric(theta_true),
    Posterior_Mean  = as.numeric(theta_post_mean),
    Bias            = as.numeric(theta_post_mean - theta_true)
  )
  
  if (!is.null(dirs)) {
    write_csv(out, file.path(dirs$tables, "dune_dirichlet_posterior_summary.csv"))
    
    # nice scatter
    if (make_plots) {
      cols <- pal_set(3, "Dark 3")
      png(file.path(dirs$figs, "dune_posterior_mean_vs_true.png"),
          width = 1200, height = 850, res = 160)
      op <- par(mar = c(5,5,4,2) + 0.1)
      plot(out$True_Proportion, out$Posterior_Mean,
           pch = 16, cex = 0.9, col = grDevices::adjustcolor(cols[1], 0.55),
           xlab = "True proportion (empirical)",
           ylab = "Posterior mean (Dirichlet)",
           main = sprintf("Dune: Dirichlet posterior mean vs true (alpha0=%.2f)", alpha0))
      abline(0, 1, lty = 2, lwd = 2, col = cols[2])
      par(op)
      dev.off()
    }
    
    # Top-12 LaTeX table
    ord <- order(counts, decreasing = TRUE)
    top <- out[ord[1:min(12, nrow(out))], ]
    
    latex_lines <- c(
      "\\begin{tabular}{lrrr}",
      "\\toprule",
      "Species & True\\_Proportion & Posterior\\_Mean & Bias \\\\",
      "\\midrule"
    )
    for (i in seq_len(nrow(top))) {
      latex_lines <- c(latex_lines, paste0(
        top$Species[i], " & ",
        fmt_num(top$True_Proportion[i], 6), " & ",
        fmt_num(top$Posterior_Mean[i], 6), " & ",
        fmt_num(top$Bias[i], 6), " \\\\"
      ))
    }
    latex_lines <- c(latex_lines, "\\bottomrule", "\\end{tabular}")
    write_lines(latex_lines, file.path(dirs$latex, "table_dune_top_species.tex"))
  }
  
  invisible(out)
}

############################
# 6) MAIN PIPELINE
############################
main <- function(CONFIG) {
  set.seed(CONFIG$seed)
  
  run_root <- file.path(CONFIG$out_base, paste0("run_", ts_string()))
  ensure_dir(run_root)
  
  dirs <- list(
    root  = run_root,
    data  = file.path(run_root, "data"),
    figs  = file.path(run_root, "figs"),
    tables= file.path(run_root, "tables"),
    latex = file.path(run_root, "latex"),
    logs  = file.path(run_root, "logs")
  )
  lapply(dirs, ensure_dir)
  
  # session info
  write_lines(capture.output(sessionInfo()), file.path(dirs$logs, "sessionInfo.txt"))
  
  # A) Pólya MC
  message("[A] Classical Pólya Monte Carlo...")
  mc <- run_polya_mc(
    a = CONFIG$polya$a, b = CONFIG$polya$b, c = CONFIG$polya$c,
    n_steps = CONFIG$polya$n_steps, R = CONFIG$polya$R, seed = CONFIG$seed
  )
  save_polya_outputs(mc, CONFIG$polya, dirs)
  
  # B) MovieLens sequential
  message("[B] MovieLens sequential (timestamp) evaluation...")
  if (!is.null(CONFIG$movielens$local_dir)) {
    u_data_path <- file.path(CONFIG$movielens$local_dir, "u.data")
    if (!file.exists(u_data_path)) stop("u.data not found at: ", u_data_path)
  } else {
    ml_store <- file.path(dirs$data, "movielens")
    ensure_dir(ml_store)
    u_data_path <- download_movielens_100k(CONFIG$movielens$zip_url, ml_store)
  }
  run_movielens_seq_polya(
    u_data_path = u_data_path,
    train_frac  = CONFIG$movielens$train_frac,
    alpha0      = CONFIG$movielens$alpha0,
    make_plots  = CONFIG$movielens$make_plots,
    dirs = dirs
  )
  
  # C) crabs
  message("[C] MASS::crabs outputs...")
  run_crabs_outputs(dirs)
  
  # D) dune Dirichlet posterior
  message("[D] vegan::dune Dirichlet posterior outputs...")
  run_dune_dirichlet_posterior(
    alpha0 = CONFIG$dune$alpha0,
    make_plots = CONFIG$dune$make_plots,
    dirs = dirs
  )
  
  # Sync latest (stable paths)
  latest_root <- file.path(CONFIG$out_base, "latest")
  message("[E] Syncing stable folder: ", latest_root)
  sync_latest(run_root, latest_root)
  
  # Final note
  write_lines(c(
    "RUN COMPLETE",
    paste0("run_root: ", run_root),
    paste0("latest_root: ", latest_root),
    "Key LaTeX tables (stable paths):",
    "  latest/latex/table_mc_bias_se_revblue.tex",
    "  latest/latex/table_movielens_sequential_eval.tex",
    "  latest/latex/table_dune_top_species.tex",
    "Key figures (stable paths):",
    "  latest/figs/polya_final_pA_hist_beta_overlay.png",
    "  latest/figs/movielens_running_mean_logpred.png",
    "  latest/figs/movielens_top_movies_observed_expected.png",
    "  latest/figs/movielens_observed_vs_expected_scatter.png",
    "  latest/figs/dune_posterior_mean_vs_true.png",
    "  latest/figs/crabs_sp_observed_vs_polya_demo.png"
  ), file.path(dirs$logs, "RUN_COMPLETE.txt"))
  
  message("DONE. Outputs saved to: ", run_root)
  message("Stable outputs available at: ", latest_root)
  
  invisible(list(run_root = run_root, latest_root = latest_root))
}

# RUN
main(CONFIG)

############################################################
# MONTE CARLO SE COMPUTATION (append after your simulation code)
############################################################

# Number of Monte Carlo replications
R_mc <- 500   # You can increase to 1000 for smoother SE

# Storage matrix:
# Each row = one replication
# Each column = final posterior mean for color j
posterior_means_mc <- matrix(NA, nrow = R_mc, ncol = 3)

for (r in 1:R_mc) {
  
  # ---- Re-run ONE full Bayesian Polya simulation ----
  
  # Reinitialize counts (same true proportions as your main simulation)
  theta_true <- c(0.3, 0.5, 0.2)
  beta <- 1
  lambda <- 2
  mu <- rep(1/3, 3)
  t_max <- 10000
  
  # Initialize counts proportional to truth
  counts <- round(theta_true * 100)   # small initial seed
  total_counts <- sum(counts)
  
  for (t in 1:t_max) {
    probs <- counts / sum(counts)
    draw <- sample(1:3, size = 1, prob = probs)
    counts[draw] <- counts[draw] + beta
  }
  
  # Posterior parameters
  alpha_post <- lambda * mu + counts
  
  # Posterior mean
  posterior_mean <- alpha_post / sum(alpha_post)
  
  # Store final posterior mean for this replication
  posterior_means_mc[r, ] <- posterior_mean
}

# ---------------------------------------------------------
# Monte Carlo mean (this is what appears in your table)
# ---------------------------------------------------------
mc_mean <- colMeans(posterior_means_mc)

# ---------------------------------------------------------
# Monte Carlo SD across replications
# ---------------------------------------------------------
mc_sd <- apply(posterior_means_mc, 2, sd)

# ---------------------------------------------------------
# Monte Carlo SE of the Monte Carlo mean
# ---------------------------------------------------------
mc_se <- mc_sd / sqrt(R_mc)

# ---------------------------------------------------------
# Bias (Monte Carlo mean - true value)
# ---------------------------------------------------------
bias <- mc_mean - theta_true

# ---------------------------------------------------------
# Final summary table for LaTeX
# ---------------------------------------------------------
mc_results <- data.frame(
  Color = 1:3,
  True = theta_true,
  Mean = round(mc_mean, 4),
  Bias = round(bias, 4),
  SE = round(mc_se, 5)
)

print(mc_results)

