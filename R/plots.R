# R/plots.R ──────────────────────────────────────────────────────────────────
# All plotly chart functions for the ISI-T2D Shiny dashboard.
# Each function returns a plotly object and handles empty / missing data.

suppressPackageStartupMessages({
  library(plotly)
  library(dplyr)
  library(tidyr)
  library(here)
})

source(here::here("R/constants.R"), local = FALSE)

# ── Shared helpers ────────────────────────────────────────────────────────────

# Standard plotly config applied to every chart
.cfg <- function(p) {
  p %>%
    plotly::layout(
      autosize = TRUE
    ) %>%
    plotly::config(
      responsive             = TRUE,
      scrollZoom             = TRUE,
      displaylogo            = FALSE,
      modeBarButtonsToRemove = list("lasso2d", "select2d", "autoScale2d"),
      toImageButtonOptions   = list(format = "svg", width = 1400, height = 900)
    )
}

# Shared layout defaults
.base_layout <- function(p, height = NULL, left_margin = 60) {
  args <- list(
    paper_bgcolor = "#fafafa",
    plot_bgcolor  = "#ffffff",
    font          = list(family = "system-ui, -apple-system, sans-serif", size = 11),
    hoverlabel    = list(bgcolor = "white", font = list(size = 12))
  )
  if (!is.null(height)) args$height <- height
  args$margin <- list(l = left_margin, r = 30, t = 40, b = 80)
  do.call(plotly::layout, c(list(p), args))
}

# Placeholder chart with a centred message
.empty_plot <- function(msg = "No data to display") {
  plotly::plot_ly() %>%
    plotly::layout(
      xaxis = list(visible = FALSE),
      yaxis = list(visible = FALSE),
      paper_bgcolor = "#fafafa",
      plot_bgcolor  = "#fafafa",
      annotations   = list(list(
        text      = msg, showarrow = FALSE,
        xref      = "paper", yref = "paper",
        x = 0.5, y = 0.5,
        font      = list(size = 16, color = "#999")
      ))
    )
}

# Compute variant row order (by gene name then variant ID)
.variant_order <- function(df) {
  df %>%
    dplyr::distinct(gene_snp_ra, nearest_gene, index_variant) %>%
    dplyr::arrange(nearest_gene, index_variant) %>%
    dplyr::pull(gene_snp_ra)
}

# Auto-scale z axis based on 99th-percentile of |beta|
.z_limit <- function(beta_vec, floor = 0.005, cap = 0.06) {
  max(floor, min(quantile(abs(beta_vec), 0.99, na.rm = TRUE) * 1.1, cap))
}

# Cell height for heatmap rows
.cell_h <- function(n) if (n <= 60) 16L else if (n <= 300) 12L else 8L

# ─────────────────────────────────────────────────────────────────────────────
# 1.  MAIN HEATMAP
# ─────────────────────────────────────────────────────────────────────────────

#' Interactive heatmap: variants (y) × IS indices (x), coloured by beta.
#' Significant cells (p < 0.05) get a black border via a second trace.
#'
#' @param df       Tidy data frame with columns: gene_snp_ra, is_index,
#'                 index_group, beta, pvalue, nearest_gene, index_variant,
#'                 t_2_d_risk_allele.
#' @param min_h    Minimum plot height in pixels (auto-expanded for many rows).
plot_heatmap <- function(df, min_h = 550) {
  if (is.null(df) || nrow(df) == 0)
    return(.empty_plot("No variants match the current filters") %>% .cfg())

  var_ord <- .variant_order(df)
  n_var   <- length(var_ord)
  ph      <- max(min_h, n_var * .cell_h(n_var) + 160L)
  zlim    <- .z_limit(df$beta)
  lm      <- min(260L, max(150L, max(nchar(var_ord)) * 6L))

  df <- df %>%
    dplyr::mutate(
      gene_snp_ra = factor(gene_snp_ra, levels = var_ord),
      is_index    = factor(is_index,    levels = IS_INDEX_ORDER),
      sig_label   = dplyr::if_else(
        pvalue < 0.001,
        " ★★",
        dplyr::if_else(pvalue < 0.05, " ★", "")
      )
    )

  plotly::plot_ly(
    data       = df,
    x          = ~is_index,
    y          = ~gene_snp_ra,
    z          = ~beta,
    type       = "heatmap",
    colorscale = BETA_COLORSCALE,
    zmin       = -zlim,
    zmax       =  zlim,
    colorbar   = list(
      title      = "β",
      thickness  = 14,
      len        = 0.55,
      tickformat = ".4f",
      tickvals   = c(-zlim, -zlim / 2, 0, zlim / 2, zlim)
    ),
    text = ~paste0(
      "<b>", gene_snp_ra, "</b>", sig_label, "<br>",
      "IS index:    ", is_index,                        "<br>",
      "Group:       ", index_group,                     "<br>",
      "β =          ", sprintf("%.5f", beta),           "<br>",
      "p =          ", formatC(pvalue, format = "e", digits = 2)
    ),
    hoverinfo = "text"
  ) %>%
    plotly::layout(
      autosize = TRUE,

      xaxis = list(
        title         = "",
        tickangle     = -48,
        tickfont      = list(size = 9),
        categoryorder = "array",
        categoryarray = as.character(IS_INDEX_ORDER),
        automargin    = TRUE,
        fixedrange    = FALSE
      ),

      yaxis = list(
        title      = "",
        tickfont   = list(size = 8),
        autorange  = "reversed",
        automargin = TRUE,
        fixedrange = FALSE
      ),

      height   = ph,
      margin   = list(l = lm, b = 120, r = 40, t = 35),
      dragmode = "pan",

      paper_bgcolor = "#fafafa",
      plot_bgcolor  = "#ffffff"
    ) %>%
    .cfg()
}

# ─────────────────────────────────────────────────────────────────────────────
# 2.  BMI ATTENUATION / DIFFERENCE HEATMAP
# ─────────────────────────────────────────────────────────────────────────────

#' Heatmap of (β_no_bmi − β_bmi): positive = BMI attenuates, negative = amplifies.
plot_bmi_diff_heatmap <- function(df_joined, min_h = 550) {
  if (is.null(df_joined) || nrow(df_joined) == 0)
    return(.empty_plot("No overlapping variants across both analyses") %>% .cfg())

  df <- df_joined %>%
    dplyr::mutate(
      beta_diff  = beta_no_bmi - beta_bmi,
      pct_change = (beta_no_bmi - beta_bmi) / (abs(beta_bmi) + 1e-10) * 100
    )

  var_ord <- .variant_order(dplyr::rename(df, beta = beta_bmi))
  n_var   <- length(var_ord)
  ph      <- max(min_h, n_var * .cell_h(n_var) + 160L)
  zlim    <- max(0.002, quantile(abs(df$beta_diff), 0.99, na.rm = TRUE) * 1.1)
  lm      <- min(260L, max(150L, max(nchar(var_ord)) * 6L))

  df <- df %>%
    dplyr::mutate(
      gene_snp_ra = factor(gene_snp_ra, levels = var_ord),
      is_index    = factor(is_index,    levels = IS_INDEX_ORDER)
    )

  plotly::plot_ly(
    data       = df,
    x          = ~is_index,
    y          = ~gene_snp_ra,
    z          = ~beta_diff,
    type       = "heatmap",
    colorscale = ATTENUATION_COLORSCALE,
    zmin       = -zlim,
    zmax       =  zlim,
    colorbar   = list(
      title      = "β diff<br>(no-BMI − BMI)",
      thickness  = 14,
      len        = 0.55,
      tickformat = ".4f"
    ),
    text = ~paste0(
      "<b>", gene_snp_ra, "</b><br>",
      "IS index:        ", is_index,                     "<br>",
      "β (BMI-adj):     ", sprintf("%.5f", beta_bmi),    "<br>",
      "β (no-BMI):      ", sprintf("%.5f", beta_no_bmi), "<br>",
      "Δβ (no-BMI−BMI): ", sprintf("%.5f", beta_diff),   "<br>",
      "% change:        ", sprintf("%.1f%%", pct_change)
    ),
    hoverinfo = "text"
  ) %>%
    plotly::layout(
      autosize = TRUE,

      xaxis = list(
        title         = "",
        tickangle     = -48,
        tickfont      = list(size = 9),
        categoryorder = "array",
        categoryarray = as.character(IS_INDEX_ORDER),
        automargin    = TRUE,
        fixedrange    = FALSE
      ),

      yaxis = list(
        title      = "",
        tickfont   = list(size = 8),
        autorange  = "reversed",
        automargin = TRUE,
        fixedrange = FALSE
      ),

      height   = ph,
      margin   = list(l = lm, b = 120, r = 40, t = 35),
      dragmode = "pan",

      paper_bgcolor = "#fafafa",
      plot_bgcolor  = "#ffffff"
    ) %>%
    .cfg()
}
# ─────────────────────────────────────────────────────────────────────────────
# 3.  VOLCANO PLOT  (single index or 3×N facet grid)
# ─────────────────────────────────────────────────────────────────────────────

plot_volcano <- function(df, selected_index = "All Indices", height = 600) {
  if (is.null(df) || nrow(df) == 0) return(.empty_plot())

  df <- df %>%
    dplyr::mutate(
      log10p = -log10(pvalue + 1e-300),
      sig    = dplyr::case_when(
        pvalue < 0.001 ~ "p < 0.001",
        pvalue < 0.05  ~ "p < 0.05",
        TRUE           ~ "NS"
      )
    )

  if (selected_index != "All Indices") {
    sub <- dplyr::filter(df, is_index == selected_index)
    if (nrow(sub) == 0) return(.empty_plot(paste("No data for:", selected_index)))

    beta_rng <- range(sub$beta, na.rm = TRUE)

    plotly::plot_ly(
      data     = sub,
      x        = ~beta,
      y        = ~log10p,
      type     = "scatter",
      mode     = "markers",
      color    = ~sig,
      colors   = SIG_COLORS,
      marker   = list(size = 8, opacity = 0.8, line = list(width = 0.5, color = "white")),
      text     = ~paste0(
        "<b>", gene_snp_ra, "</b><br>",
        "β = ", sprintf("%.5f", beta),  "<br>",
        "p = ", formatC(pvalue, format = "e", digits = 2)
      ),
      hoverinfo = "text"
    ) %>%
      plotly::layout(
        title  = list(text = selected_index, font = list(size = 13)),
        xaxis  = list(title = "Effect size (β)", zeroline = TRUE,
                      zerolinecolor = "#bbb", zerolinewidth = 1.5),
        yaxis  = list(title = "-log₁₀(p-value)"),
        height = height,
        shapes = list(
          list(type = "line", xref = "paper", x0 = 0, x1 = 1,
               y0 = -log10(0.05), y1 = -log10(0.05),
               line = list(dash = "dash", color = "#fc8d59", width = 1.2)),
          list(type = "line", xref = "paper", x0 = 0, x1 = 1,
               y0 = -log10(0.001), y1 = -log10(0.001),
               line = list(dash = "dot",  color = "#d73027", width = 1.2))
        ),
        legend = list(title = list(text = "Significance")),
        paper_bgcolor = "#fafafa",
        plot_bgcolor  = "#ffffff"
      ) %>%
      .cfg()

  } else {
    # Faceted 3-column grid
    indices <- intersect(IS_INDEX_ORDER, unique(df$is_index))
    n_cols  <- 3L
    n_rows  <- ceiling(length(indices) / n_cols)

    plots <- lapply(seq_along(indices), function(i) {
      idx <- indices[[i]]
      sub <- dplyr::filter(df, is_index == idx)
      is_left <- ((i - 1L) %% n_cols) == 0L
      is_first <- i == 1L

      plotly::plot_ly(
        data      = sub,
        x         = ~beta,
        y         = ~log10p,
        type      = "scatter",
        mode      = "markers",
        color     = ~sig,
        colors    = SIG_COLORS,
        showlegend = is_first,
        marker    = list(size = 5, opacity = 0.75),
        text      = ~paste0("<b>", gene_snp_ra, "</b><br>β=",
                            sprintf("%.4f", beta), "<br>p=",
                            formatC(pvalue, format = "e", digits = 1)),
        hoverinfo = "text"
      ) %>%
        plotly::layout(
          annotations = list(list(
            text      = idx, showarrow = FALSE,
            xref      = "paper", yref = "paper",
            x = 0.5, y = 1.05, font = list(size = 9, color = "#333")
          )),
          xaxis = list(title = if (i > (n_rows - 1L) * n_cols) "β" else "",
                       zeroline = TRUE, zerolinecolor = "#ccc", tickfont = list(size = 8)),
          yaxis = list(title = if (is_left) "-log₁₀p" else "",
                       tickfont = list(size = 8))
        )
    })

    plotly::subplot(plots, nrows = n_rows, shareX = FALSE, shareY = FALSE,
                    titleX = TRUE, titleY = TRUE, margin = 0.05) %>%
      plotly::layout(
        height        = max(600, n_rows * 210),
        paper_bgcolor = "#fafafa",
        plot_bgcolor  = "#ffffff",
        legend        = list(title = list(text = "Sig."), x = 1.01, y = 1,
                             font = list(size = 10))
      ) %>%
      .cfg()
  }
}

# ─────────────────────────────────────────────────────────────────────────────
# 4.  EFFECT PROFILE  (forest-style dot plot for one variant)
# ─────────────────────────────────────────────────────────────────────────────

#' Show β for all 21 IS indices for a selected variant.
#' Both BMI-adjusted and unadjusted are plotted as separate symbols.
plot_forest <- function(df_bmi, df_no_bmi, variant, height = 560) {
  if (is.null(variant) || nchar(trimws(variant)) == 0)
    return(.empty_plot("Select a variant in the dropdown above"))

  d_bmi    <- dplyr::filter(df_bmi,    gene_snp_ra == variant)
  d_no_bmi <- dplyr::filter(df_no_bmi, gene_snp_ra == variant)

  if (nrow(d_bmi) + nrow(d_no_bmi) == 0)
    return(.empty_plot(paste("Variant not found:", variant)))

  d_bmi$analysis    <- "BMI-Adjusted"
  d_no_bmi$analysis <- "BMI-Unadjusted"

  df_all <- dplyr::bind_rows(d_bmi, d_no_bmi) %>%
    dplyr::mutate(
      is_index = factor(is_index, levels = rev(IS_INDEX_ORDER)),
      # index group shapes for y-axis region bands
      group_num = match(index_group, c("Fasting", "OGTT,0-120", "OGTT,0-30-120"))
    )

  cols <- c("BMI-Adjusted" = "#2c7bb6", "BMI-Unadjusted" = "#d7191c")

  # Error bars from SE column if present; otherwise ±0.5 * |beta| (cosmetic)
  has_se <- "se" %in% names(df_all)

  p <- plotly::plot_ly(
    data     = df_all,
    x        = ~beta,
    y        = ~is_index,
    color    = ~analysis,
    colors   = cols,
    type     = "scatter",
    mode     = "markers",
    marker   = list(size = 11, opacity = 0.85,
                    line = list(width = 1.5, color = "white")),
    symbol   = ~analysis,
    symbols  = c("circle", "diamond"),
    text     = ~paste0(
      "<b>", is_index, "</b>  (", index_group, ")<br>",
      "Analysis: ", analysis, "<br>",
      "β = ",       sprintf("%.5f", beta),  "<br>",
      "p = ",       formatC(pvalue, format = "e", digits = 2),
      if_else(pvalue < 0.05, "<br><b>★ p < 0.05</b>", "")
    ),
    hoverinfo = "text"
  )

  if (has_se) {
    for (ana in c("BMI-Adjusted", "BMI-Unadjusted")) {
      sub <- dplyr::filter(df_all, analysis == ana)
      p <- p %>% plotly::add_trace(
        data      = sub,
        x         = ~beta,
        y         = ~is_index,
        error_x   = list(array = sub$se, color = cols[[ana]], thickness = 1.5, width = 4),
        type      = "scatter",
        mode      = "markers",
        marker    = list(opacity = 0),
        showlegend = FALSE,
        hoverinfo = "skip",
        inherit   = FALSE
      )
    }
  }

  p %>%
    plotly::layout(
      title  = list(text = paste("Effect profile:", variant), font = list(size = 13)),
      xaxis  = list(
        title        = "Effect size (β)",
        zeroline     = TRUE,
        zerolinecolor = "#999",
        zerolinewidth = 2
      ),
      yaxis  = list(title = "", tickfont = list(size = 10)),
      height = height,
      legend = list(title = list(text = "Analysis"), x = 1.01, y = 0.5),
      margin = list(l = 150, r = 130, t = 50, b = 60),
      paper_bgcolor = "#fafafa",
      plot_bgcolor  = "#ffffff",
      shapes = list(
        list(type = "line", x0 = 0, x1 = 0,
             y0 = 0, y1 = 1, yref = "paper",
             line = list(color = "#aaa", width = 1.5, dash = "dot"))
      )
    ) %>%
    .cfg()
}

# ─────────────────────────────────────────────────────────────────────────────
# 5.  β SCATTER: BMI-adjusted vs BMI-unadjusted
# ─────────────────────────────────────────────────────────────────────────────

plot_beta_scatter <- function(df_joined, height = 580) {
  if (is.null(df_joined) || nrow(df_joined) == 0) return(.empty_plot())

  df <- df_joined %>%
    dplyr::mutate(
      sig_status = dplyr::case_when(
        pval_bmi < 0.05 & pval_no_bmi < 0.05 ~ "Sig. in both",
        pval_bmi < 0.05                        ~ "Sig. BMI-adj only",
        pval_no_bmi < 0.05                     ~ "Sig. no-BMI only",
        TRUE                                   ~ "Not significant"
      )
    )

  cor_val  <- round(cor(df$beta_bmi, df$beta_no_bmi, use = "complete.obs"), 3)
  br       <- range(c(df$beta_bmi, df$beta_no_bmi), na.rm = TRUE) * 1.05

  plotly::plot_ly(
    data     = df,
    x        = ~beta_bmi,
    y        = ~beta_no_bmi,
    type     = "scatter",
    mode     = "markers",
    color    = ~index_group,
    colors   = INDEX_GROUP_COLORS,
    symbol   = ~sig_status,
    symbols  = c("circle", "square", "diamond", "cross"),
    marker   = list(size = 6, opacity = 0.72, line = list(width = 0.3, color = "white")),
    text     = ~paste0(
      "<b>", gene_snp_ra, "</b><br>",
      "IS Index:      ", is_index, "  (", index_group, ")<br>",
      "β (BMI-adj):   ", sprintf("%.5f", beta_bmi),    "<br>",
      "β (no-BMI):    ", sprintf("%.5f", beta_no_bmi), "<br>",
      "p (BMI-adj):   ", formatC(pval_bmi,    format = "e", digits = 2), "<br>",
      "p (no-BMI):    ", formatC(pval_no_bmi, format = "e", digits = 2)
    ),
    hoverinfo = "text"
  ) %>%
    # Diagonal identity line
    plotly::add_trace(
      x = br, y = br,
      type = "scatter", mode = "lines",
      line = list(color = "#888", dash = "dash", width = 1.5),
      showlegend = FALSE, hoverinfo = "none", inherit = FALSE
    ) %>%
    plotly::add_annotations(
      x = br[2] * 0.9, y = br[1] * 1.4,
      text      = paste0("r = ", cor_val),
      showarrow = FALSE,
      font      = list(size = 13, color = "#333")
    ) %>%
    plotly::layout(
      xaxis  = list(title = "β (BMI-adjusted)",   zeroline = TRUE, zerolinecolor = "#ccc"),
      yaxis  = list(title = "β (BMI-unadjusted)", zeroline = TRUE, zerolinecolor = "#ccc"),
      height = height,
      legend = list(title = list(text = "Index Group / Status"), x = 1.01),
      margin = list(l = 70, r = 160, t = 40, b = 70),
      paper_bgcolor = "#fafafa",
      plot_bgcolor  = "#ffffff"
    ) %>%
    .cfg()
}

# ─────────────────────────────────────────────────────────────────────────────
# 6.  PARALLEL COORDINATES
# ─────────────────────────────────────────────────────────────────────────────

plot_parallel <- function(df, height = 580) {
  if (is.null(df) || nrow(df) == 0) return(.empty_plot())

  df_wide <- df %>%
    dplyr::select(gene_snp_ra, nearest_gene, is_index, beta) %>%
    tidyr::pivot_wider(names_from = is_index, values_from = beta) %>%
    dplyr::mutate(gene_num = as.numeric(factor(nearest_gene)))

  idx_cols   <- intersect(IS_INDEX_ORDER, names(df_wide))
  all_betas  <- unlist(df_wide[idx_cols], use.names = FALSE)
  global_rng <- range(all_betas, na.rm = TRUE)

  dimensions <- lapply(idx_cols, function(col) {
    list(
      label  = col,
      values = df_wide[[col]],
      range  = global_rng
    )
  })

  plotly::plot_ly(
    df_wide,
    type       = "parcoords",
    line       = list(
      color      = ~gene_num,
      colorscale = "Viridis",
      showscale  = FALSE,
      opacity    = 0.45,
      reversescale = FALSE
    ),
    dimensions = dimensions
  ) %>%
    plotly::layout(
      height        = height,
      margin        = list(l = 70, r = 70, t = 60, b = 40),
      paper_bgcolor = "#fafafa"
    ) %>%
    .cfg()
}

# ─────────────────────────────────────────────────────────────────────────────
# 7.  GENE SUMMARY BAR CHART
# ─────────────────────────────────────────────────────────────────────────────

plot_gene_summary <- function(df, top_n = 30L, height = 500) {
  if (is.null(df) || nrow(df) == 0) return(.empty_plot())

  gene_tbl <- df %>%
    dplyr::group_by(nearest_gene) %>%
    dplyr::summarise(
      n_sig         = sum(pvalue < 0.05,  na.rm = TRUE),
      n_sig_strict  = sum(pvalue < 0.001, na.rm = TRUE),
      best_p        = min(pvalue,         na.rm = TRUE),
      mean_abs_beta = mean(abs(beta),     na.rm = TRUE),
      n_variants    = dplyr::n_distinct(gene_snp_ra),
      .groups = "drop"
    ) %>%
    dplyr::arrange(dplyr::desc(n_sig), best_p) %>%
    dplyr::slice_head(n = top_n) %>%
    dplyr::mutate(nearest_gene = factor(nearest_gene, levels = rev(nearest_gene)))

  plotly::plot_ly(
    data        = gene_tbl,
    x           = ~n_sig,
    y           = ~nearest_gene,
    type        = "bar",
    orientation = "h",
    marker      = list(
      color        = ~best_p,
      colorscale   = list(list(0, "#d73027"), list(0.5, "#fee090"), list(1, "#74add1")),
      reversescale = TRUE,
      showscale    = TRUE,
      colorbar     = list(title = "Min p", len = 0.5, tickformat = ".2e")
    ),
    hovertext   = ~paste0(
      "<b>", nearest_gene, "</b><br>",
      "Sig. assoc. (p<0.05):  ", n_sig,         "<br>",
      "Sig. assoc. (p<0.001): ", n_sig_strict,  "<br>",
      "Best p-value:          ", formatC(best_p, format = "e", digits = 2), "<br>",
      "# Variants:            ", n_variants,    "<br>",
      "Mean |β|:              ", sprintf("%.5f", mean_abs_beta)
    ),
    hoverinfo   = "text"
  ) %>%
    plotly::layout(
      xaxis  = list(title = "# Significant associations (p < 0.05)"),
      yaxis  = list(title = "", tickfont = list(size = 9)),
      height = height,
      margin = list(l = 120, r = 60, t = 35, b = 60),
      paper_bgcolor = "#fafafa",
      plot_bgcolor  = "#ffffff"
    ) %>%
    .cfg()
}

# ─────────────────────────────────────────────────────────────────────────────
# 8.  β DISTRIBUTION: violin + box per IS index
# ─────────────────────────────────────────────────────────────────────────────

plot_index_violin <- function(df, height = 500) {
  if (is.null(df) || nrow(df) == 0) return(.empty_plot())

  df <- df %>%
    dplyr::mutate(is_index = factor(is_index, levels = IS_INDEX_ORDER))

  plotly::plot_ly(
    data     = df,
    x        = ~is_index,
    y        = ~beta,
    color    = ~index_group,
    colors   = INDEX_GROUP_COLORS,
    type     = "violin",
    box      = list(visible = TRUE),
    meanline = list(visible = TRUE, color = "white"),
    points   = "outliers",
    pointpos = 0,
    hovertemplate = "<b>%{x}</b><br>β = %{y:.5f}<extra>%{fullData.name}</extra>"
  ) %>%
    plotly::layout(
      xaxis   = list(title = "", tickangle = -48, tickfont = list(size = 9)),
      yaxis   = list(title = "Effect size (β)", zeroline = TRUE, zerolinecolor = "#bbb"),
      height  = height,
      violinmode = "group",
      shapes  = list(
        list(type = "line", xref = "paper", x0 = 0, x1 = 1,
             y0 = 0, y1 = 0,
             line = list(color = "#bbb", dash = "dot", width = 1.5))
      ),
      legend  = list(title = list(text = "Index Group")),
      margin  = list(l = 65, r = 30, t = 35, b = 120),
      paper_bgcolor = "#fafafa",
      plot_bgcolor  = "#ffffff"
    ) %>%
    .cfg()
}

# ─────────────────────────────────────────────────────────────────────────────
# 9.  BMI ATTENUATION BAR CHART  (top variants by BMI mediation)
# ─────────────────────────────────────────────────────────────────────────────

plot_attenuation_bar <- function(df_joined, top_n = 30L, height = 500) {
  if (is.null(df_joined) || nrow(df_joined) == 0) return(.empty_plot())

  df <- df_joined %>%
    dplyr::mutate(
      pct_mediated = (beta_no_bmi - beta_bmi) / (abs(beta_no_bmi) + 1e-10) * 100
    ) %>%
    dplyr::group_by(gene_snp_ra) %>%
    dplyr::summarise(
      mean_pct  = mean(pct_mediated, na.rm = TRUE),
      max_abs   = max(abs(beta_bmi), na.rm = TRUE),
      best_p    = min(pval_bmi,      na.rm = TRUE),
      .groups   = "drop"
    ) %>%
    dplyr::arrange(dplyr::desc(abs(mean_pct))) %>%
    dplyr::slice_head(n = top_n) %>%
    dplyr::mutate(
      gene_snp_ra = factor(gene_snp_ra, levels = rev(gene_snp_ra)),
      direction   = dplyr::if_else(mean_pct > 0, "BMI attenuates", "BMI amplifies")
    )

  plotly::plot_ly(
    data        = df,
    x           = ~mean_pct,
    y           = ~gene_snp_ra,
    type        = "bar",
    orientation = "h",
    color       = ~direction,
    colors      = c("BMI attenuates" = "#762a83", "BMI amplifies" = "#1b7837"),
    hovertext   = ~paste0(
      "<b>", gene_snp_ra, "</b><br>",
      "Mean % mediated by BMI: ", sprintf("%.1f%%", mean_pct), "<br>",
      "Max |β| (BMI-adj):      ", sprintf("%.5f", max_abs),    "<br>",
      "Best p (BMI-adj):       ", formatC(best_p, format = "e", digits = 2)
    ),
    hoverinfo   = "text"
  ) %>%
    plotly::layout(
      xaxis  = list(title = "Mean % change in β: (no-BMI − BMI) / |no-BMI| × 100"),
      yaxis  = list(title = "", tickfont = list(size = 9)),
      height = height,
      shapes = list(
        list(type = "line", x0 = 0, x1 = 0,
             y0 = 0, y1 = 1, yref = "paper",
             line = list(color = "#aaa", width = 1.5))
      ),
      legend = list(title = list(text = "Direction")),
      margin = list(l = 200, r = 60, t = 35, b = 60),
      paper_bgcolor = "#fafafa",
      plot_bgcolor  = "#ffffff"
    ) %>%
    .cfg()
}
