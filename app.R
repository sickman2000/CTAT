# =============================================================================
# CRISPR Trace Shiny App
# =============================================================================

library(rstudioapi)
if (rstudioapi::isAvailable()) {
  setwd(dirname(getActiveDocumentContext()$path))
}

library(shiny)
library(sangerseqR)
library(Biostrings)
library(nnls)

source("R/tide_functions.R", local = TRUE)

# -----------------------------------------------------------------------------
ui <- fluidPage(
  titlePanel("CRISPR Trace Analysis"),
  sidebarLayout(
    sidebarPanel(
      width = 3,
      fileInput("wt_file", "Wild-type .ab1", accept = c(".ab1", ".AB1")),
      fileInput("crispr_file", "CRISPR .ab1", accept = c(".ab1", ".AB1")),
      textInput("grna", "gRNA (20 bp, without PAM)",
                value = "TTCGAGGGGAAGACATCGTT"),
      textAreaInput("master_seq", "Reference sequence (optional, DNA)",
                    value = "", rows = 3,
                    placeholder = "FASTA or plain DNA — for gRNA validation and when the WT trace is too short"),
      checkboxInput("force_partial", "Fallback mode (analyze even if the gRNA is missing from the WT)", FALSE),
      helpText("Use only if the WT barely covers the cut site; results are then indicative only."),
      numericInput("upstream", "Upstream (bases)", value = 25, min = 5, max = 50),
      numericInput("window_size", "Downstream (bases)", value = 100, min = 30, max = 200),
      numericInput("max_indel", "Max. indel size", value = 10, min = 5, max = 20),
      actionButton("run", "Start analysis", class = "btn-primary", width = "100%"),
      hr(),
      conditionalPanel(
        condition = "output.has_results",
        h4("Top indels for trace comparison"),
        uiOutput("indel_checkboxes"),
        checkboxInput("show_mixture", "Mixture of selected indels", TRUE),
        checkboxInput("show_cut", "Show cut site", TRUE)
      )
    ),
    mainPanel(
      width = 9,
      tabsetPanel(
        id = "tabs",
        tabPanel(
          "Overview",
          br(),
          verbatimTextOutput("log"),
          hr(),
          fluidRow(
            column(6, plotOutput("spectrum_plot", height = "350px")),
            column(6, tableOutput("summary_table"))
          )
        ),
        tabPanel(
          "Trace comparison",
          br(),
          plotOutput("overlay_plot", height = "320px"),
          hr(),
          h4("Indel model vs. CRISPR (downstream)"),
          plotOutput("indel_trace_plot", height = "420px"),
          helpText(
            "Blue = WT (unshifted). Red = CRISPR. ",
            "Green dashed line = cut site. ",
            "Colored lines = expected WT pattern for each indel ",
            "(downstream shifted only)."
          )
        ),
        tabPanel(
          "4-channel WT vs. indel traces",
          br(),
          h4("Synchronized four-channel and sequence comparison: wild-type vs. indel model"),
          helpText("Displays the WT and each indel model in four channels, with synchronized sequences arranged one below the other."),
          uiOutput("four_channel_indel_plots_ui")
        ),
        tabPanel(
          "Statistics",
          br(),
          downloadButton("dl_csv", "Frequencies as CSV"),
          br(), br(),
          tableOutput("freq_table")
        )
      )
    )
  )
)

# -----------------------------------------------------------------------------
server <- function(input, output, session) {

  rv <- reactiveValues(
    pipeline = NULL,
    log = "Ready. Load AB1 files and start the analysis."
  )

  append_log <- function(msg) {
    rv$log <- paste(rv$log, msg, sep = "\n")
  }

  observeEvent(input$run, {
    req(input$wt_file, input$crispr_file, input$grna)

    grna <- toupper(gsub("[^ACGTacgt]", "", input$grna))
    if (nchar(grna) < 15) {
      showNotification("gRNA is too short or invalid.", type = "error")
      return()
    }

    rv$log <- "Starting analysis ..."
    withProgress(message = "CRISPR analysis", value = 0, {

      tryCatch({
        incProgress(0.15, detail = "Load AB1 files and locate gRNA")
        master_raw <- if (!is.null(input$master_seq) && nzchar(input$master_seq))
          input$master_seq else NULL

        data <- initialize_tide_analysis(
          wt_file       = input$wt_file$datapath,
          crispr_file   = input$crispr_file$datapath,
          grna_seq      = grna,
          master_seq    = master_raw,
          force_partial = isTRUE(input$force_partial)
        )
        append_log(sprintf(
          "gRNA identity (%.1f%%, %s). Cut site: base %d | Mode: %s",
          data$identity * 100,
          if (data$is_reverse) "Antisense" else "Sense",
          data$cut_site,
          data$mode
        ))
        if (!is.null(data$master)) {
          m <- data$master
          append_log(sprintf(
            "Reference: gRNA at %d–%d (cut %d, %.1f%%) | WT covers reference %d–%d (cut covered by WT: %s)",
            m$grna_start, m$grna_end, m$cut_site, m$identity * 100,
            m$wt_cov_start, m$wt_cov_end,
            if (isTRUE(m$wt_covers_cut)) "yes" else "NO"
          ))
        }
        for (w in data$warnings) append_log(paste("WARNING:", w))

        incProgress(0.25, detail = "WT window")
        trace_data <- get_analysis_traces(
          data,
          window_size = as.integer(input$window_size),
          upstream    = as.integer(input$upstream)
        )

        incProgress(0.35, detail = "Trace alignment (peak-to-peak spline)")
        aligned <- align_crispr_traces_affine(data, trace_data)
        append_log(sprintf(
          "Peak-to-peak spline: offset=%d | velocity=%.4f | peak pairs=%d | spline=%s | pre r=%.3f | near r=%.3f | post r=%.3f",
          aligned$optimal_offset,
          aligned$warp_velocity,
          aligned$n_anchors,
          if (isTRUE(aligned$spline_used)) "yes" else "approximation",
          aligned$pre_cut_correlation,
          aligned$near_cut_correlation,
          aligned$post_cut_correlation
        ))
        if (is.finite(aligned$near_cut_correlation) &&
            aligned$near_cut_correlation < 0.70 &&
            aligned$pre_cut_correlation > 0.9) {
          append_log(sprintf(
            "Note: near-cut correlation (%.3f) is low; the peak grid may not yet be optimal.",
            aligned$near_cut_correlation
          ))
        }
        append_log(
          "Info: Peak maxima follow product length and should overlap; A/C/G/T proportions may differ downstream of the cut."
        )

        incProgress(0.15, detail = "Amplitude normalization")
        aligned <- normalize_amplitudes(trace_data, aligned)

        incProgress(0.15, detail = "Deconvolution")
        results <- deconvolve_indels(
          trace_data, aligned,
          max_indel = as.integer(input$max_indel)
        )
        append_log(sprintf(
          "Editing efficiency: %.1f%% | WT: %.1f%% | R2: %.3f | bp_step: %.3f data points/base",
          results$efficiency,
          results$frequencies["Indel_0"],
          results$r_squared,
          results$bp_step
        ))

        rv$pipeline <- list(
          tide_data  = data,
          trace_data = trace_data,
          aligned    = aligned,
          results    = results
        )
        showNotification("Analysis complete.", type = "message")
      }, error = function(e) {
        append_log(paste("ERROR:", e$message))
        showNotification(e$message, type = "error", duration = 10)
        rv$pipeline <- NULL
      })
    })
  })

  output$has_results <- reactive(!is.null(rv$pipeline))
  outputOptions(output, "has_results", suspendWhenHidden = FALSE)

  output$log <- renderText(rv$log)

  # ---- Indel checkboxes ----
  output$indel_checkboxes <- renderUI({
    req(rv$pipeline)
    freq <- rv$pipeline$results$frequencies
    top <- sort(freq, decreasing = TRUE)
    top <- top[top >= 1]
    if (length(top) == 0) top <- sort(freq, decreasing = TRUE)[1:5]

    choices <- names(top)
    selected <- names(top)[seq_len(min(4, length(top)))]

    checkboxGroupInput(
      "selected_indels",
      label = NULL,
      choices = setNames(
        choices,
        sprintf("%s (%.1f%%)", choices, unname(top[choices]))
      ),
      selected = selected
    )
  })

  # ---- Summary ----
  output$summary_table <- renderTable({
    req(rv$pipeline)
    r <- rv$pipeline$results
    a <- rv$pipeline$aligned
    d <- rv$pipeline$tide_data
    data.frame(
      Parameter = c(
        "Cut site (WT)", "Direction", "gRNA identity", "Mode",
        "Pre-cut correlation", "Near-cut correlation", "Post-cut correlation",
        "Warp velocity", "Peak anchors", "Offset",
        "bp_step (dp/Base)", "Editing efficiency", "Wild-type", "R²"
      ),
      Value = c(
        d$cut_site,
        if (d$is_reverse) "Antisense" else "Sense",
        sprintf("%.1f%%", d$identity * 100),
        if (!is.null(d$mode)) d$mode else "normal",
        sprintf("%.3f", a$pre_cut_correlation),
        sprintf("%.3f", a$near_cut_correlation),
        sprintf("%.3f", a$post_cut_correlation),
        sprintf("%.4f", a$warp_velocity),
        a$n_anchors,
        a$optimal_offset,
        sprintf("%.2f", r$bp_step),
        sprintf("%.1f%%", r$efficiency),
        sprintf("%.1f%%", r$frequencies["Indel_0"]),
        sprintf("%.3f", r$r_squared)
      )
    )
  }, striped = TRUE, spacing = "s")

  # ---- Spectrum barplot ----
  output$spectrum_plot <- renderPlot({
    req(rv$pipeline)
    freq <- rv$pipeline$results$frequencies
    cols <- ifelse(names(freq) == "Indel_0", "#2c7bb6", "#d7191c")
    barplot(
      freq,
      las = 2, col = cols, border = NA,
      main = sprintf(
        "Indel spectrum (efficiency %.1f%% | R2=%.2f)",
        rv$pipeline$results$efficiency,
        rv$pipeline$results$r_squared
      ),
      ylab = "Fraction (%)",
      cex.names = 0.75
    )
    abline(h = 0)
  })

  # ---- Frequency table + download ----
  output$freq_table <- renderTable({
    req(rv$pipeline)
    freq <- rv$pipeline$results$frequencies
    data.frame(
      Indel = names(freq),
      Fraction_Percent = round(unname(freq), 2),
      row.names = NULL
    )
  }, striped = TRUE)

  output$dl_csv <- downloadHandler(
    filename = function() paste0("crispr_indels_", Sys.Date(), ".csv"),
    content = function(file) {
      req(rv$pipeline)
      freq <- rv$pipeline$results$frequencies
      write.csv(
        data.frame(Indel = names(freq), Fraction_Percent = unname(freq)),
        file, row.names = FALSE
      )
    }
  )

  # ---- Overlay: WT vs. CRISPR full window ----
  output$overlay_plot <- renderPlot({
    req(rv$pipeline)
    td <- rv$pipeline$trace_data
    al <- rv$pipeline$aligned
    cut_pos <- td$cut_trace_index - td$start_trace_index + 1L

    wt_sig     <- rowSums(td$wt_window_trace)
    crispr_sig <- rowSums(al$crispr_window_trace)

    ylim <- range(c(wt_sig, crispr_sig), finite = TRUE)

    plot(wt_sig, type = "l", col = "#2c7bb6", lwd = 1.8, ylim = ylim,
         xlab = "Data point within window", ylab = "Intensity (sum of A+C+G+T)",
         main = sprintf(
           "WT (blue) vs CRISPR (red) | Pre r=%.3f | Near r=%.3f | Post r=%.3f",
           al$pre_cut_correlation,
           al$near_cut_correlation,
           al$post_cut_correlation
         ))
    lines(crispr_sig, col = "#d7191c", lwd = 1.4)
    if (isTRUE(input$show_cut)) {
      abline(v = cut_pos, col = "#1a9641", lty = 2, lwd = 2)
    }
    legend("topright",
           c("WT", "CRISPR", "Cut site"),
           col = c("#2c7bb6", "#d7191c", "#1a9641"),
           lwd = 2, lty = c(1, 1, 2), bty = "n")
  })

  # ---- Indel model traces vs CRISPR ----
  output$indel_trace_plot <- renderPlot({
    req(rv$pipeline, input$selected_indels)

    res <- rv$pipeline$results
    cut_pos <- res$cut_pos
    n <- nrow(res$wt_mat)

    wt_sum     <- rowSums(res$wt_mat)
    crispr_sum <- rowSums(res$crispr_mat)

    indel_names <- input$selected_indels
    n_ind <- length(indel_names)
    pal <- if (n_ind <= 8) {
      c("#e41a1c", "#984ea3", "#ff7f00", "#a65628",
        "#f781bf", "#999999", "#66c2a5", "#fc8d62")[seq_len(n_ind)]
    } else {
      grDevices::hcl.colors(n_ind, "Dark 3")
    }
    names(pal) <- indel_names

    ylim <- range(c(wt_sum, crispr_sum), finite = TRUE)

    for (nm in indel_names) {
      if (!is.null(res$model_matrices[[nm]])) {
        ylim <- range(c(ylim, rowSums(res$model_matrices[[nm]])), finite = TRUE)
      }
    }

    plot(NA, xlim = c(1, n), ylim = ylim,
         xlab = "Data point within window", ylab = "Intensity",
         main = "Indel models (downstream shift) vs. CRISPR")

    lines(wt_sum, col = "#2c7bb6", lwd = 2)
    lines(crispr_sum, col = "#d7191c", lwd = 1.6)

    for (nm in indel_names) {
      mat <- res$model_matrices[[nm]]
      if (is.null(mat)) next
      lines(rowSums(mat), col = pal[nm], lwd = 1.5, lty = 1)
    }

    if (isTRUE(input$show_mixture) && length(indel_names) >= 1) {
      mix <- matrix(0, nrow = n, ncol = 4)
      total_w <- sum(res$frequencies[indel_names])
      if (total_w > 0) {
        for (nm in indel_names) {
          w <- res$frequencies[nm] / total_w
          mix <- mix + w * res$model_matrices[[nm]]
        }
        lines(rowSums(mix), col = "black", lwd = 2, lty = 3)
      }
    }

    if (isTRUE(input$show_cut)) {
      abline(v = cut_pos, col = "#1a9641", lty = 2, lwd = 2)
    }

    leg <- c("WT", "CRISPR", indel_names)
    leg_col <- c("#2c7bb6", "#d7191c", unname(pal[indel_names]))
    leg_lty <- rep(1, length(leg))
    if (isTRUE(input$show_mixture) && length(indel_names) >= 1) {
      leg <- c(leg, "Mixture")
      leg_col <- c(leg_col, "black")
      leg_lty <- c(leg_lty, 3)
    }
    if (isTRUE(input$show_cut)) {
      leg <- c(leg, "Cut site")
      leg_col <- c(leg_col, "#1a9641")
      leg_lty <- c(leg_lty, 2)
    }

    legend("topright", legend = leg, col = leg_col, lty = leg_lty,
           lwd = 2, bty = "n", cex = 0.85)
  })

  # ---- Dynamic UI: actual WT/CRISPR traces, plus indel model and sequence alignment ----
  output$four_channel_indel_plots_ui <- renderUI({
    req(rv$pipeline, input$selected_indels)
    indel_names <- input$selected_indels
    if (!length(indel_names)) return(p("Please select at least one indel on the left."))

    do.call(tagList, lapply(seq_along(indel_names), function(i) {
      tagList(
        wellPanel(
          h4(sprintf("Indel model: %s (%.2f%%)", indel_names[i],
                     rv$pipeline$results$frequencies[indel_names[i]])),
          plotOutput(paste0("four_channel_comp_plot_", i), height = "1100px"),
          h5("WT and CRISPR sequences: alignment at the target region"),
          verbatimTextOutput(paste0("seq_text_comparison_", i))
        ),
        br()
      )
    }))
  })

observe({
  req(rv$pipeline, input$selected_indels)
  indel_names <- input$selected_indels
  pipe <- rv$pipeline
  lapply(seq_along(indel_names), function(i) {
    local({
      idx <- i
      nm <- indel_names[idx]
      plotname <- paste0("four_channel_comp_plot_", idx)
      textname <- paste0("seq_text_comparison_", idx)

      output[[plotname]] <- renderPlot({
        req(rv$pipeline)
        p <- rv$pipeline
        res <- p$results
        wt_mat <- res$wt_mat
        crispr_mat <- res$crispr_mat
        model_mat <- res$model_matrices[[nm]]

        validate(need(!is.null(crispr_mat), "No aligned CRISPR trace is available."))
        validate(need(!is.null(model_mat), paste("No model data for", nm)))

        channel_cols <- c(A = "#009E73", C = "#0072B2",
                          G = "#222222", T = "#D55E00")
        cut_pos <- max(1, min(res$cut_pos, nrow(wt_mat)))

        plot_one <- function(mat, title_text, show_xlab = TRUE) {
          ymax <- max(mat, na.rm = TRUE)
          if (!is.finite(ymax) || ymax <= 0) ymax <- 1
          plot(
            NA,
            xlim = c(1, nrow(mat)),
            ylim = c(0, ymax * 1.18),
            xlab = if (show_xlab) "Relative trace data point (mapped to WT window)" else "",
            ylab = "Signal intensity",
            main = title_text
          )
          for (ch in c("A", "C", "G", "T")) {
            lines(seq_len(nrow(mat)), mat[, ch], col = channel_cols[ch], lwd = 1.1)
          }
          if (isTRUE(input$show_cut)) {
            abline(v = cut_pos, col = "#1a9641", lty = 2, lwd = 1.5)
          }
          legend(
            "topright",
            legend = c("A", "C", "G", "T",
                       if (isTRUE(input$show_cut)) "Cut site"),
            col = c(unname(channel_cols),
                    if (isTRUE(input$show_cut)) "#1a9641"),
            lty = c(rep(1, 4), if (isTRUE(input$show_cut)) 2),
            lwd = c(rep(1.5, 4), if (isTRUE(input$show_cut)) 1.5),
            bty = "n", horiz = TRUE, cex = 0.8
          )
        }

        oldpar <- par(no.readonly = TRUE)
        on.exit(par(oldpar))

        # Mehr Abstand zwischen den Panels + Platz für den äußeren Titel
        par(
          mfrow = c(3, 1),
          mar   = c(4.5, 4.5, 4.0, 1.2),   # bottom, left, top, right – hier die Abstände steuern
          oma   = c(0, 0, 2.8, 0),         # oberer Außenrand für den Gesamttitel
          mgp   = c(2.3, 0.7, 0)
        )

        plot_one(wt_mat, "WT – measured 4-channel trace", show_xlab = FALSE)
        plot_one(
          crispr_mat,
          sprintf(
            "CRISPR – peak warping (pre r=%.3f | near r=%.3f | post r=%.3f)",
            p$aligned$pre_cut_correlation,
            p$aligned$near_cut_correlation,
            p$aligned$post_cut_correlation
          ),
          show_xlab = FALSE
        )
        plot_one(
          model_mat,
          sprintf("%s – WT model trace with %s bp shift from cut site",
                  nm, sub("Indel_", "", nm)),
          show_xlab = TRUE
        )

        mtext(
          sprintf("%s | fraction in NNLS model: %.2f%%",
                  nm, res$frequencies[nm]),
          side = 3, outer = TRUE, line = 0.9, cex = 1.05, font = 2
        )
      }, res = 110)

      output[[textname]] <- renderText({
        req(rv$pipeline)
        p <- rv$pipeline
        td <- p$trace_data
        wt_full <- toupper(as.character(p$tide_data$wt_seq))
        start_b <- max(1L, as.integer(td$start_base))
        end_b <- min(nchar(wt_full), as.integer(td$end_base))
        wt_window <- substr(wt_full, start_b, end_b)
        wt_chars <- strsplit(wt_window, "", fixed = TRUE)[[1]]
        n_wt <- length(wt_chars)
        indel_shift <- as.integer(sub("^Indel_", "", nm))
        # Convention: the cut lies after the WT base at cut_site.
        # This alignment is a synthetic model, not an actual CRISPR read.
        cut_rel <- as.integer(p$tide_data$cut_site) - start_b + 1L
        cut_rel <- max(0L, min(cut_rel, n_wt))
        prefix <- if (cut_rel > 0L) wt_chars[seq_len(cut_rel)] else character(0)
        suffix_start <- cut_rel + 1L
        if (indel_shift < 0L) {
          # Deletion: removed WT bases are represented by gaps in the CRISPR alignment.
          n_del <- min(abs(indel_shift), max(0L, n_wt - cut_rel))
          del_idx <- if (n_del > 0L) seq.int(cut_rel + 1L, cut_rel + n_del) else integer(0)
          suffix_idx <- cut_rel + n_del + 1L
          deleted <- if (n_del > 0L) wt_chars[del_idx] else character(0)
          suffix <- if (suffix_idx <= n_wt) wt_chars[seq.int(suffix_idx, n_wt)] else character(0)
          awt_chars <- c(prefix, deleted, suffix)
          acr_chars <- c(prefix, rep("-", n_del), suffix)
          shift_desc <- sprintf("Model: %d bp deletion", indel_shift)
        } else if (indel_shift > 0L) {
          # A pure shift model does not know the inserted bases: N = unknown.
          suffix <- if (suffix_start <= n_wt) wt_chars[seq.int(suffix_start, n_wt)] else character(0)
          awt_chars <- c(prefix, rep("-", indel_shift), suffix)
          acr_chars <- c(prefix, rep("N", indel_shift), suffix)
          shift_desc <- sprintf("Model: +%d bp insertion (inserted bases unknown)", indel_shift)
        } else {
          awt_chars <- wt_chars
          acr_chars <- wt_chars
          shift_desc <- "Model: no indel (WT)"
        }
        awt <- paste0(awt_chars, collapse = "")
        acr <- paste0(acr_chars, collapse = "")
        matchline <- ifelse(awt_chars == acr_chars, "|",
                            ifelse(awt_chars == "-" | acr_chars == "-", " ", "."))
        matchline <- paste0(matchline, collapse = "")
        wrap_alignment <- function(x, width = 80L) {
          if (!nchar(x)) return(character(0))
          starts <- seq.int(1L, nchar(x), by = width)
          vapply(starts, function(s) substr(x, s, min(nchar(x), s + width - 1L)),
                 character(1))
        }
        lines_wt <- wrap_alignment(awt)
        lines_cr <- wrap_alignment(acr)
        lines_m <- wrap_alignment(matchline)
        chunks <- unlist(lapply(seq_along(lines_wt), function(j) c(
          sprintf("WT     %s", lines_wt[j]),
          sprintf("       %s", lines_m[j]),
          sprintf("CRISPR %s", lines_cr[j]), ""
        )), use.names = FALSE)
        paste(c(
          sprintf("%s | estimated NNLS fraction: %.2f%%",
                  nm, p$results$frequencies[nm]),
          shift_desc,
          "Artificial model alignment: | = match, . = difference, gap = deletion/insertion",
          sprintf("WT window: bases %d–%d; model cut after WT base %d.",
                  start_b, end_b, as.integer(p$tide_data$cut_site)),
          "The alignment and model trace are calculated; they are not actual CRISPR sequencing data.",
          if (indel_shift > 0L) "N represents the unknown inserted bases." else NULL,
          if (indel_shift < 0L && abs(indel_shift) > max(0L, n_wt - cut_rel))
            "WARNING: The window does not contain enough downstream bases for the full modeled deletion." else NULL,
          "",
          chunks
        ), collapse = "\n")
      })
    })
  })
})
}

# -----------------------------------------------------------------------------
shinyApp(ui, server)
