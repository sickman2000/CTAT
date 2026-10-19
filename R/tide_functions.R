# =============================================================================
# CRISPR Trace Analysis – Function Library (sequential peak-to-peak spline)
# =============================================================================
# Dependencies: sangerseqR, Biostrings, nnls
# =============================================================================
#
# Core principle:
#   Peak position = product length. Capillary separation varies locally.
#   Pair successive peaks by index after a common start, fit a spline
#   CRISPR_xpos ~ f(WT_xpos), and warp the CRISPR trace onto the WT grid.
#
# Alignment strategy (v6):
# 1. Map one common start base (pre-cut) via sequence alignment.
# 2. Pair WT peak[start+k] with CRISPR peak[start_cr+k] for the whole window.
# 3. smooth.spline over all pairs (captures local mobility differences).
# 4. Dense warp of every WT data point through the spline.
# 5. Residual integer offset on pre-cut only.
# 6. Indel models shift intensity patterns by N peak intervals on this grid.
# =============================================================================

#' @title Initialize CRISPR Trace analysis
#' @description Load AB1 files, locate gRNA, compute SpCas9 cut site.
#'   Optional master_seq: reference DNA used to validate gRNA presence and
#'   to recover the cut site when the WT trace is truncated.
#' @param force_partial If TRUE and gRNA is missing from WT, try to place
#'   the cut via master-WT alignment (last-resort mode).
initialize_tide_analysis <- function(wt_file,
                                     crispr_file,
                                     grna_seq,
                                     min_identity = 0.8,
                                     master_seq = NULL,
                                     force_partial = FALSE) {

  wt_sanger     <- sangerseqR::readsangerseq(wt_file)
  crispr_sanger <- sangerseqR::readsangerseq(crispr_file)

  wt_string     <- sangerseqR::primarySeq(wt_sanger)
  crispr_string <- sangerseqR::primarySeq(crispr_sanger)

  grna    <- Biostrings::DNAString(grna_seq)
  grna_rc <- Biostrings::reverseComplement(grna)

  warnings_list <- character(0)
  master_info   <- NULL
  mode          <- "normal"

  # ----- Optional master sequence -----
  if (!is.null(master_seq) && nzchar(master_seq)) {
    master_clean <- toupper(gsub("[^ACGTacgt]", "", master_seq))
    if (nchar(master_clean) < nchar(grna_seq)) {
      stop("Master sequence is shorter than the gRNA.")
    }
    master_dna <- Biostrings::DNAString(master_clean)

    # gRNA must be in master (or RC)
    m_aln    <- Biostrings::pairwiseAlignment(grna, master_dna, type = "local")
    m_aln_rc <- Biostrings::pairwiseAlignment(grna_rc, master_dna, type = "local")
    m_id     <- Biostrings::nmatch(m_aln)    / Biostrings::nchar(grna)
    m_id_rc  <- Biostrings::nmatch(m_aln_rc) / Biostrings::nchar(grna)

    master_reverse <- FALSE
    if (m_id_rc > m_id) {
      m_aln <- m_aln_rc
      m_id  <- m_id_rc
      master_reverse <- TRUE
    }
    if (m_id < min_identity) {
      stop(sprintf(
        "gRNA not found in the master sequence (best identity %.1f%%). Please check the gRNA and master sequence.",
        m_id * 100
      ))
    }

    m_range <- Biostrings::subject(m_aln)
    grna_start_master <- BiocGenerics::start(m_range)
    grna_end_master   <- BiocGenerics::end(m_range)
    if (master_reverse) {
      cut_master <- grna_start_master + 3L
    } else {
      cut_master <- grna_end_master - 3L
    }

    # Align WT primary sequence to master to see coverage
    wt_vs_m <- Biostrings::pairwiseAlignment(wt_string, master_dna, type = "global-local")
    wt_cov_start <- BiocGenerics::start(Biostrings::subject(wt_vs_m))
    wt_cov_end   <- BiocGenerics::end(Biostrings::subject(wt_vs_m))
    wt_covers_cut <- (cut_master >= wt_cov_start && cut_master <= wt_cov_end)
    wt_covers_grna <- (grna_start_master >= wt_cov_start &&
                       grna_end_master   <= wt_cov_end)

    # Map master cut position → WT base index via the alignment
    cut_in_wt_via_master <- NA_integer_
    wt_aln_str <- as.character(Biostrings::alignedPattern(wt_vs_m))
    m_aln_str  <- as.character(Biostrings::alignedSubject(wt_vs_m))
    wt_pos <- BiocGenerics::start(Biostrings::pattern(wt_vs_m)) - 1L
    m_pos  <- BiocGenerics::start(Biostrings::subject(wt_vs_m)) - 1L
    chars_wt <- strsplit(wt_aln_str, "")[[1]]
    chars_m  <- strsplit(m_aln_str, "")[[1]]
    for (i in seq_along(chars_m)) {
      if (chars_wt[i] != "-") wt_pos <- wt_pos + 1L
      if (chars_m[i]  != "-") m_pos  <- m_pos + 1L
      if (m_pos == cut_master && chars_wt[i] != "-") {
        cut_in_wt_via_master <- wt_pos
        break
      }
    }

    master_info <- list(
      sequence          = master_clean,
      grna_start        = as.integer(grna_start_master),
      grna_end          = as.integer(grna_end_master),
      cut_site          = as.integer(cut_master),
      is_reverse        = master_reverse,
      identity          = m_id,
      wt_cov_start      = as.integer(wt_cov_start),
      wt_cov_end        = as.integer(wt_cov_end),
      wt_covers_cut     = wt_covers_cut,
      wt_covers_grna    = wt_covers_grna,
      cut_in_wt         = cut_in_wt_via_master
    )
  }

  # ----- Locate gRNA in WT trace -----
  alignment    <- Biostrings::pairwiseAlignment(grna, wt_string, type = "local")
  alignment_rc <- Biostrings::pairwiseAlignment(grna_rc, wt_string, type = "local")

  identity_score    <- Biostrings::nmatch(alignment)    / Biostrings::nchar(grna)
  identity_score_rc <- Biostrings::nmatch(alignment_rc) / Biostrings::nchar(grna)

  is_reverse <- FALSE
  if (identity_score_rc > identity_score) {
    alignment      <- alignment_rc
    identity_score <- identity_score_rc
    is_reverse     <- TRUE
  }

  if (identity_score >= min_identity) {
    wt_range         <- Biostrings::subject(alignment)
    grna_start_in_wt <- BiocGenerics::start(wt_range)
    grna_end_in_wt   <- BiocGenerics::end(wt_range)
    if (is_reverse) {
      cut_site_wt <- grna_start_in_wt + 3L
    } else {
      cut_site_wt <- grna_end_in_wt - 3L
    }
    mode <- "normal"
  } else {
    # gRNA not found in WT – try master-based recovery
    if (!is.null(master_info) && isTRUE(master_info$wt_covers_cut) &&
        !is.na(master_info$cut_in_wt)) {
      cut_site_wt      <- master_info$cut_in_wt
      grna_start_in_wt <- max(1L, cut_site_wt - 17L)
      grna_end_in_wt   <- min(length(wt_string), cut_site_wt + 3L)
      is_reverse       <- master_info$is_reverse
      identity_score   <- master_info$identity
      mode <- "master_rescue"
      warnings_list <- c(warnings_list, sprintf(
        "gRNA not found in the WT trace (%.1f%%). Cut site mapped to WT base %d using the master sequence.",
        max(identity_score, identity_score_rc) * 100, cut_site_wt
      ))
    } else if (isTRUE(force_partial) && !is.null(master_info) &&
               !is.na(master_info$cut_in_wt)) {
      # Desperate mode: WT may not fully cover cut; use mapped position anyway
      cut_site_wt <- master_info$cut_in_wt
      n_wt <- length(wt_string)
      if (cut_site_wt < 1L || cut_site_wt > n_wt) {
        # Cut outside WT – clamp to end and shrink window later
        cut_site_wt <- max(1L, min(n_wt, cut_site_wt))
        warnings_list <- c(warnings_list,
          "EMERGENCY: Cut site is outside or at the edge of the WT trace; results are highly uncertain.")
      }
      grna_start_in_wt <- max(1L, cut_site_wt - 17L)
      grna_end_in_wt   <- min(n_wt, cut_site_wt + 3L)
      is_reverse       <- master_info$is_reverse
      identity_score   <- max(identity_score, identity_score_rc)
      mode <- "force_partial"
      warnings_list <- c(warnings_list, sprintf(
        "EMERGENCY MODE: gRNA missing from the WT trace (%.1f%%). Analysis uses master-sequence mapping (cut ~ base %d); results are indicative only.",
        identity_score * 100, cut_site_wt
      ))
    } else if (!is.null(master_info) && !isTRUE(master_info$wt_covers_cut)) {
      stop(sprintf(
        paste0(
          "gRNA not found in the WT trace (%.1f%%). The master sequence indicates that the WT covers positions %d–%d, ",
          "the cut site is at master position %d, so the WT trace is too short. ",
          "Load a longer WT trace or enable 'Emergency Mode'."
        ),
        max(identity_score, identity_score_rc) * 100,
        master_info$wt_cov_start, master_info$wt_cov_end, master_info$cut_site
      ))
    } else {
      stop(sprintf(
        "gRNA not found. Best identity: %.1f%%%s",
        max(identity_score, identity_score_rc) * 100,
        if (is.null(master_info)) " (optionally upload a master sequence for diagnosis)." else "."
      ))
    }
  }

  # Sanity: cut inside WT
  n_wt <- length(wt_string)
  if (cut_site_wt < 5L || cut_site_wt > n_wt - 5L) {
    warnings_list <- c(warnings_list, sprintf(
      "Cut site (base %d) is near the edge of the WT trace (%d bases); the downstream window is limited.",
      cut_site_wt, n_wt
    ))
  }

  list(
    wt_sanger     = wt_sanger,
    crispr_sanger = crispr_sanger,
    wt_seq        = wt_string,
    crispr_seq    = crispr_string,
    cut_site      = as.integer(cut_site_wt),
    grna_start    = as.integer(grna_start_in_wt),
    grna_end      = as.integer(grna_end_in_wt),
    is_reverse    = is_reverse,
    identity      = identity_score,
    master        = master_info,
    mode          = mode,
    warnings      = warnings_list
  )
}

#' @title Extract WT analysis window (upstream + downstream of cut)
get_analysis_traces <- function(tide_data, window_size = 100L, upstream = 25L) {

  wt_obj   <- tide_data$wt_sanger
  cut_base <- tide_data$cut_site

  peak_pos_matrix   <- sangerseqR::peakPosMatrix(wt_obj)
  wt_peak_positions <- rowMeans(peak_pos_matrix, na.rm = TRUE)

  n_bases <- length(sangerseqR::primarySeq(wt_obj))
  if (length(wt_peak_positions) > n_bases) {
    wt_peak_positions <- wt_peak_positions[seq_len(n_bases)]
  }

  start_base <- max(1L, cut_base - upstream)
  end_base   <- min(cut_base + window_size, length(wt_peak_positions))

  # Cas9 cuts BETWEEN nucleotides. cut_base is the last base before the cut
  # (e.g. position 17 of the spacer). Place the cut marker at the midpoint
  # between peak[cut_base] and peak[cut_base+1], not on a peak maximum.
  peak_before <- wt_peak_positions[cut_base]
  if (cut_base < length(wt_peak_positions)) {
    peak_after <- wt_peak_positions[cut_base + 1L]
    cut_trace_index <- as.integer(round((peak_before + peak_after) / 2))
  } else {
    cut_trace_index <- as.integer(round(peak_before))
  }
  start_trace_idx <- as.integer(round(wt_peak_positions[start_base]))
  end_trace_index <- as.integer(round(wt_peak_positions[end_base]))

  n_trace_rows <- nrow(sangerseqR::traceMatrix(wt_obj))
  start_trace_idx <- max(1L, min(start_trace_idx, n_trace_rows))
  end_trace_index <- max(start_trace_idx, min(end_trace_index, n_trace_rows))
  cut_trace_index <- max(start_trace_idx, min(cut_trace_index, n_trace_rows))

  wt_traces_window <- sangerseqR::traceMatrix(wt_obj)[
    start_trace_idx:end_trace_index, , drop = FALSE
  ]
  colnames(wt_traces_window) <- c("A", "C", "G", "T")

  list(
    cut_trace_index   = cut_trace_index,
    start_trace_index = start_trace_idx,
    end_trace_index   = end_trace_index,
    start_base        = start_base,
    end_base          = end_base,
    wt_window_trace   = wt_traces_window,
    peak_positions    = wt_peak_positions
  )
}

#' @title Map one or more WT base positions to the CRISPR read
#' @return Integer vector of same length as wt_bases (NA where unmapped)
map_bases_wt_to_crispr <- function(tide_data, wt_bases) {

  wt_seq     <- sangerseqR::primarySeq(tide_data$wt_sanger)
  crispr_seq <- sangerseqR::primarySeq(tide_data$crispr_sanger)

  aln <- Biostrings::pairwiseAlignment(wt_seq, crispr_seq, type = "global-local")

  wt_str     <- as.character(Biostrings::alignedPattern(aln))
  crispr_str <- as.character(Biostrings::alignedSubject(aln))

  wt_pos     <- BiocGenerics::start(Biostrings::pattern(aln)) - 1L
  crispr_pos <- BiocGenerics::start(Biostrings::subject(aln)) - 1L

  chars_wt     <- strsplit(wt_str, "")[[1]]
  chars_crispr <- strsplit(crispr_str, "")[[1]]

  max_wt <- length(wt_seq)
  mapping <- rep(NA_integer_, max_wt)

  for (i in seq_along(chars_wt)) {
    if (chars_wt[i] != "-")     wt_pos     <- wt_pos + 1L
    if (chars_crispr[i] != "-") crispr_pos <- crispr_pos + 1L
    if (chars_wt[i] != "-" && chars_crispr[i] != "-" &&
        wt_pos >= 1L && wt_pos <= max_wt) {
      mapping[wt_pos] <- crispr_pos
    }
  }

  as.integer(mapping[wt_bases])
}

#' @title Map a single WT base (compatibility wrapper)
map_base_wt_to_crispr <- function(tide_data, wt_base) {
  map_bases_wt_to_crispr(tide_data, as.integer(wt_base))[1]
}

#' @title Align CRISPR traces via sequential peak-to-peak spline
#'
#' Capillary separation varies locally over the run. Method:
#'   1. Map one common start base (pre-cut) WT -> CRISPR via sequence.
#'   2. Pair peaks by sequential index from that start over the WHOLE window:
#'      WT peak[start+k] <-> CRISPR peak[start_cr+k]  (pre- and post-cut).
#'      Successive peaks = successive product lengths; no sequence identity needed.
#'   3. Spline regression: CRISPR_xpos ~ f(WT_xpos) on all pairs.
#'   4. Apply the spline as dense warp to every WT data point.
#'   5. Small residual integer offset on pre-cut only.
align_crispr_traces_affine <- function(tide_data,
                                       wt_trace_info,
                                       search_padding = 40L) {

  crispr_obj  <- tide_data$crispr_sanger
  cut_base_wt <- tide_data$cut_site

  crispr_peaks   <- rowMeans(sangerseqR::peakPosMatrix(crispr_obj), na.rm = TRUE)
  n_crispr_rows  <- nrow(sangerseqR::traceMatrix(crispr_obj))
  n_crispr_bases <- length(sangerseqR::primarySeq(crispr_obj))
  if (length(crispr_peaks) > n_crispr_bases) {
    crispr_peaks <- crispr_peaks[seq_len(n_crispr_bases)]
  }

  start_base_wt <- wt_trace_info$start_base
  end_base_wt   <- wt_trace_info$end_base
  wt_peaks_all  <- wt_trace_info$peak_positions

  window_length <- nrow(wt_trace_info$wt_window_trace)
  cut_pos <- wt_trace_info$cut_trace_index - wt_trace_info$start_trace_index + 1L
  cut_pos <- max(2L, min(cut_pos, window_length - 1L))

  # ----- 1. Common start (sequence mapping, pre-cut, reliable) -----
  start_cr <- map_base_wt_to_crispr(tide_data, start_base_wt)
  if (is.na(start_cr) || start_cr < 1L || start_cr > length(crispr_peaks)) {
    stop("Start-base mapping failed.")
  }

  # ----- 2. Sequential peak pairs over the whole window -----
  n_win  <- end_base_wt - start_base_wt + 1L
  max_k  <- n_win - 1L
  max_k_cr <- length(crispr_peaks) - start_cr
  max_k <- min(max_k, max_k_cr)
  if (max_k < 5L) stop("Too few sequential peak pairs.")

  ks     <- 0:max_k
  wt_idx <- start_base_wt + ks
  cr_idx <- start_cr + ks

  ok <- wt_idx >= 1L & wt_idx <= length(wt_peaks_all) &
        cr_idx >= 1L & cr_idx <= length(crispr_peaks)

  fit_wt <- as.numeric(wt_peaks_all[wt_idx[ok]])
  fit_cr <- as.numeric(crispr_peaks[cr_idx[ok]])

  ord <- order(fit_wt)
  fit_wt <- fit_wt[ord]
  fit_cr <- fit_cr[ord]
  keep <- c(TRUE, diff(fit_wt) > 0.5)
  fit_wt <- fit_wt[keep]
  fit_cr <- fit_cr[keep]

  for (k in seq_along(fit_cr)[-1]) {
    if (fit_cr[k] <= fit_cr[k - 1]) fit_cr[k] <- fit_cr[k - 1] + 0.5
  }

  n_fit <- length(fit_wt)
  if (n_fit < 6L) stop("Too few peak pairs remain after filtering.")

  i0 <- max(1L, n_fit - 3L)
  vel <- (fit_cr[n_fit] - fit_cr[i0]) / max(1e-6, fit_wt[n_fit] - fit_wt[i0])
  if (!is.finite(vel) || vel < 0.4 || vel > 2.5) {
    vel <- stats::median(diff(fit_cr) / pmax(1e-6, diff(fit_wt)))
    if (!is.finite(vel) || vel < 0.4 || vel > 2.5) vel <- 1.0
  }

  # ----- 3. Spline on ALL sequential pairs (local mobility correction) -----
  spline_fit <- tryCatch(
    stats::smooth.spline(x = fit_wt, y = fit_cr, spar = 0.35, all.knots = TRUE),
    error = function(e) NULL
  )
  use_spline <- !is.null(spline_fit)

  # ----- 4. Dense warp -----
  start_abs <- wt_trace_info$start_trace_index
  end_abs   <- wt_trace_info$end_trace_index
  wt_dp_abs <- as.numeric(seq.int(start_abs, end_abs))

  if (use_spline) {
    warp_cr <- as.numeric(stats::predict(spline_fit, wt_dp_abs)$y)
    below <- wt_dp_abs < fit_wt[1]
    above <- wt_dp_abs > fit_wt[n_fit]
    if (any(below)) {
      v0 <- (fit_cr[min(2L, n_fit)] - fit_cr[1]) /
            max(1e-6, fit_wt[min(2L, n_fit)] - fit_wt[1])
      warp_cr[below] <- fit_cr[1] + v0 * (wt_dp_abs[below] - fit_wt[1])
    }
    if (any(above)) {
      warp_cr[above] <- fit_cr[n_fit] + vel * (wt_dp_abs[above] - fit_wt[n_fit])
    }
  } else {
    warp_cr <- stats::approx(fit_wt, fit_cr, xout = wt_dp_abs, rule = 2)$y
  }

  # ----- 5. Residual integer offset (pre-cut) -----
  wt_pre <- rowSums(wt_trace_info$wt_window_trace[seq_len(cut_pos), , drop = FALSE])
  crispr_trace_full <- sangerseqR::traceMatrix(crispr_obj)

  best_cor    <- -Inf
  best_offset <- 0L
  best_mat    <- NULL

  for (off in -search_padding:search_padding) {
    sample_pos <- warp_cr + off
    if (min(sample_pos) < 1 || max(sample_pos) > n_crispr_rows) next

    resampled <- matrix(0, nrow = window_length, ncol = 4L)
    for (ch in 1:4) {
      resampled[, ch] <- stats::approx(
        x = seq_len(n_crispr_rows),
        y = crispr_trace_full[, ch],
        xout = sample_pos,
        rule = 2
      )$y
    }

    test_pre <- rowSums(resampled[seq_len(cut_pos), , drop = FALSE])
    if (stats::sd(test_pre) < 1e-6) next

    this_cor <- stats::cor(wt_pre, test_pre)
    if (this_cor > best_cor) {
      best_cor    <- this_cor
      best_offset <- off
      best_mat    <- resampled
    }
  }

  if (is.null(best_mat)) stop("Alignment failed: no valid offset found.")
  colnames(best_mat) <- c("A", "C", "G", "T")

  post_start <- min(cut_pos + 1L, window_length)
  if (post_start < window_length) {
    wt_post     <- rowSums(wt_trace_info$wt_window_trace[post_start:window_length, , drop = FALSE])
    crispr_post <- rowSums(best_mat[post_start:window_length, , drop = FALSE])
    post_cor <- if (stats::sd(wt_post) > 1e-6 && stats::sd(crispr_post) > 1e-6)
      stats::cor(wt_post, crispr_post) else NA_real_
  } else post_cor <- NA_real_

  near_end <- min(window_length,
                  cut_pos + max(20L, as.integer(0.15 * (window_length - cut_pos))))
  if (near_end > cut_pos + 5L) {
    wt_near     <- rowSums(wt_trace_info$wt_window_trace[(cut_pos + 1L):near_end, , drop = FALSE])
    crispr_near <- rowSums(best_mat[(cut_pos + 1L):near_end, , drop = FALSE])
    near_cor <- if (stats::sd(wt_near) > 1e-6 && stats::sd(crispr_near) > 1e-6)
      stats::cor(wt_near, crispr_near) else NA_real_
  } else near_cor <- NA_real_

  cut_base_crispr <- map_base_wt_to_crispr(tide_data, cut_base_wt)

  list(
    optimal_offset        = best_offset,
    scale_factor          = vel,
    scale_estimate        = vel,
    pre_cut_correlation   = best_cor,
    post_cut_correlation  = post_cor,
    near_cut_correlation  = near_cor,
    crispr_window_trace   = best_mat,
    cut_base_crispr       = cut_base_crispr,
    warp_velocity         = vel,
    n_anchors             = n_fit,
    n_grid_peaks          = n_fit,
    spline_used           = use_spline
  )
}

#' @title Normalize CRISPR amplitudes using pre-cut region
normalize_amplitudes <- function(wt_trace_info, crispr_aligned_data) {

  cut_pos <- wt_trace_info$cut_trace_index - wt_trace_info$start_trace_index + 1L
  cut_pos <- min(
    cut_pos,
    nrow(wt_trace_info$wt_window_trace),
    nrow(crispr_aligned_data$crispr_window_trace)
  )

  wt_pre_mean     <- mean(wt_trace_info$wt_window_trace[seq_len(cut_pos), ])
  crispr_pre_mean <- mean(crispr_aligned_data$crispr_window_trace[seq_len(cut_pos), ])

  if (!is.finite(crispr_pre_mean) || crispr_pre_mean < 1e-8) {
    amp_scale <- 1
  } else {
    amp_scale <- wt_pre_mean / crispr_pre_mean
  }
  crispr_aligned_data$crispr_window_trace <-
    crispr_aligned_data$crispr_window_trace * amp_scale
  crispr_aligned_data$amp_scale <- amp_scale
  crispr_aligned_data
}

#' @title Shift matrix downstream of cut site using peak-position distances
#'
#' When peak_positions and start_base/cut_base are supplied, the shift in
#' data-points is taken from the actual distance between the cut-base peak
#' and the peak of (cut_base + shift_bp). Falls back to bp_step * shift_bp.
shift_matrix_downstream <- function(mat, shift_bp, bp_step, cut_pos,
                                    peak_positions = NULL,
                                    start_base = NULL,
                                    cut_base = NULL,
                                    start_trace_index = NULL) {
  n <- nrow(mat)
  shifted <- mat

  if (shift_bp == 0L) return(shifted)

  # Prefer actual peak distance when available
  shift_dp <- NA_integer_
  if (!is.null(peak_positions) && !is.null(cut_base) &&
      !is.null(start_base) && !is.null(start_trace_index)) {
    target_base <- cut_base + shift_bp
    if (target_base >= 1L && target_base <= length(peak_positions) &&
        cut_base >= 1L && cut_base <= length(peak_positions)) {
      # data-point distance between the two peaks, relative to the window
      dp_cut    <- peak_positions[cut_base]
      dp_target <- peak_positions[target_base]
      shift_dp  <- as.integer(round(dp_target - dp_cut))
    }
  }
  if (!is.finite(shift_dp)) {
    shift_dp <- as.integer(round(shift_bp * bp_step))
  }

  if (shift_dp == 0L) return(shifted)

  downstream_len <- n - cut_pos + 1L

  if (shift_dp > 0L) {
    shifted[cut_pos:n, ] <- 0
    if (shift_dp < downstream_len) {
      shifted[(cut_pos + shift_dp):n, ] <- mat[cut_pos:(n - shift_dp), ]
    }
  } else {
    abs_shift <- abs(shift_dp)
    shifted[cut_pos:n, ] <- 0
    if (abs_shift < downstream_len) {
      shifted[cut_pos:(n - abs_shift), ] <- mat[(cut_pos + abs_shift):n, ]
    }
  }
  shifted
}

#' @title Robust peak-spacing estimate for the analysis window
estimate_bp_step <- function(wt_trace_info) {
  # Robust inter-peak spacing (data points per base) from WT AB1 peaks.
  # Uses median of pairwise diffs after MAD-based outlier rejection.
  # Prefers the downstream half of the window (most relevant for indel shifts).

  peaks <- wt_trace_info$peak_positions
  start_b <- as.integer(wt_trace_info$start_base)
  end_b   <- as.integer(wt_trace_info$end_base)

  idx <- start_b:end_b
  idx <- idx[idx >= 1L & idx <= length(peaks)]
  if (length(idx) < 3L) {
    # Fallback from window length in data points / number of bases
    n_dp <- nrow(wt_trace_info$wt_window_trace)
    n_b  <- max(1L, end_b - start_b)
    return(max(5, min(25, n_dp / n_b)))
  }

  p <- as.numeric(peaks[idx])
  d <- diff(p)
  d <- d[is.finite(d) & d > 0]

  # MAD outlier filter (keep values within median +/- 3*MAD)
  robust_filter <- function(x) {
    if (length(x) < 3L) return(x)
    med <- stats::median(x)
    mad <- stats::mad(x, constant = 1)
    if (!is.finite(mad) || mad < 1e-6) {
      # no variation – keep all positive finite
      return(x[x > 0])
    }
    lo <- med - 3 * mad
    hi <- med + 3 * mad
    x[x >= max(1, lo) & x <= hi]
  }

  d <- robust_filter(d)
  if (length(d) < 2L) {
    n_dp <- nrow(wt_trace_info$wt_window_trace)
    n_b  <- max(1L, end_b - start_b)
    return(max(5, min(25, n_dp / n_b)))
  }

  # Prefer downstream half (after ~30 % of the window bases)
  mid <- start_b + as.integer((end_b - start_b) * 0.30)
  down_idx <- idx[idx >= mid]
  if (length(down_idx) >= 5L) {
    d_down <- robust_filter(diff(as.numeric(peaks[down_idx])))
    d_down <- d_down[is.finite(d_down) & d_down > 0]
    if (length(d_down) >= 3L) {
      step <- as.numeric(stats::median(d_down))
      if (is.finite(step) && step >= 4 && step <= 30) return(step)
    }
  }

  step <- as.numeric(stats::median(d))
  if (!is.finite(step) || step < 4 || step > 30) {
    n_dp <- nrow(wt_trace_info$wt_window_trace)
    n_b  <- max(1L, end_b - start_b)
    step <- n_dp / n_b
  }
  max(4, min(30, step))
}

#' @title NNLS deconvolution of indel spectrum
deconvolve_indels <- function(wt_trace_info, crispr_aligned_data, max_indel = 10L) {

  wt_mat     <- wt_trace_info$wt_window_trace
  crispr_mat <- crispr_aligned_data$crispr_window_trace

  if (nrow(wt_mat) != nrow(crispr_mat)) {
    stop("Window lengths do not match.")
  }

  n_dp <- nrow(wt_mat)
  y    <- as.vector(crispr_mat)

  cut_pos <- wt_trace_info$cut_trace_index - wt_trace_info$start_trace_index + 1L
  cut_pos <- max(2L, min(cut_pos, n_dp - 1L))

  bp_step <- estimate_bp_step(wt_trace_info)
  if (!is.finite(bp_step) || bp_step < 4 || bp_step > 30) {
    n_b <- max(1L, wt_trace_info$end_base - wt_trace_info$start_base)
    bp_step <- max(4, min(30, n_dp / n_b))
  }

  # cut_base: base whose peak is closest to the cut data-point
  peaks <- wt_trace_info$peak_positions
  cut_base_seq <- which.min(abs(peaks - wt_trace_info$cut_trace_index))
  if (length(cut_base_seq) == 0L || is.na(cut_base_seq[1])) {
    cut_base_seq <- wt_trace_info$start_base +
      as.integer(round((cut_pos - 1) / max(bp_step, 1)))
  }
  cut_base_seq <- as.integer(cut_base_seq[1])
  # Clamp into the analysis window
  cut_base_seq <- max(wt_trace_info$start_base,
                      min(wt_trace_info$end_base, cut_base_seq))

  indels <- -max_indel:max_indel
  X <- matrix(0, nrow = length(y), ncol = length(indels))
  colnames(X) <- paste0("Indel_", indels)

  model_matrices <- vector("list", length(indels))
  names(model_matrices) <- colnames(X)

  for (i in seq_along(indels)) {
    shifted <- shift_matrix_downstream(
      wt_mat, indels[i], bp_step, cut_pos,
      peak_positions     = peaks,
      start_base         = wt_trace_info$start_base,
      cut_base           = as.integer(cut_base_seq),
      start_trace_index  = wt_trace_info$start_trace_index
    )
    X[, i] <- as.vector(shifted)
    model_matrices[[i]] <- shifted
  }

  fit <- nnls::nnls(X, y)

  coefs <- fit$x
  names(coefs) <- colnames(X)
  total <- sum(coefs)
  if (total < 1e-8) stop("Deconvolution failed.")

  percentages <- (coefs / total) * 100
  wt_pct      <- unname(percentages["Indel_0"])
  efficiency  <- 100 - wt_pct
  r_squared   <- 1 - sum(fit$residuals^2) / sum((y - mean(y))^2)

  list(
    frequencies     = percentages,
    efficiency      = as.numeric(efficiency),
    r_squared       = r_squared,
    coefficients    = coefs,
    residuals       = fit$residuals,
    cut_pos         = cut_pos,
    bp_step         = bp_step,
    model_matrices  = model_matrices,
    wt_mat          = wt_mat,
    crispr_mat      = crispr_mat
  )
}

#' @title Full pipeline
run_tide_pipeline <- function(wt_file,
                              crispr_file,
                              grna_seq,
                              window_size = 100L,
                              upstream = 25L,
                              max_indel = 10L,
                              min_identity = 0.8,
                              master_seq = NULL,
                              force_partial = FALSE) {

  data <- initialize_tide_analysis(
    wt_file, crispr_file, grna_seq, min_identity,
    master_seq = master_seq, force_partial = force_partial
  )
  trace_data <- get_analysis_traces(data, window_size, upstream)
  aligned <- align_crispr_traces_affine(data, trace_data)
  aligned <- normalize_amplitudes(trace_data, aligned)
  results <- deconvolve_indels(trace_data, aligned, max_indel)

  list(
    tide_data   = data,
    trace_data  = trace_data,
    aligned     = aligned,
    results     = results
  )
}

#' @title Build reconstructed mixture signal from top indels
reconstruct_mixture <- function(results, top_n = 5L) {

  freq <- results$frequencies
  top  <- sort(freq, decreasing = TRUE)
  top  <- top[top > 0.5]
  if (length(top) > top_n) top <- top[seq_len(top_n)]

  n <- nrow(results$wt_mat)
  mixture <- matrix(0, nrow = n, ncol = 4L)
  colnames(mixture) <- c("A", "C", "G", "T")

  for (nm in names(top)) {
    w <- unname(top[nm]) / 100
    mixture <- mixture + w * results$model_matrices[[nm]]
  }

  list(mixture = mixture, components = top)
}

library(ggplot2)
library(tidyr)
library(dplyr)

plot_indel_trace_4channel <- function(chromatogram_data, start_pos, end_pos, indel_label = "") {
  sub_data <- chromatogram_data %>%
    filter(Position >= start_pos & Position <= end_pos)

  plot_data <- sub_data %>%
    pivot_longer(cols = c(A, C, G, T), names_to = "Channel", values_to = "Intensity")

  channel_colors <- c("A" = "green", "C" = "blue", "G" = "black", "T" = "red")

  p <- ggplot(plot_data, aes(x = Position, y = Intensity, color = Channel)) +
    geom_line(linewidth = 0.8) +
    scale_color_manual(values = channel_colors) +
    theme_minimal() +
    labs(
      title = paste("4-channel trace for indel:", indel_label),
      subtitle = paste("Positions:", start_pos, "to", end_pos),
      x = "Position",
      y = "Signal intensity",
      color = "Channel"
    ) +
    theme(
      panel.grid.minor = element_blank(),
      legend.position = "bottom"
    )

  if ("Base" %in% colnames(sub_data)) {
    p <- p + geom_text(
      data = sub_data,
      aes(x = Position, y = max(plot_data$Intensity, na.rm = TRUE) * 1.05, label = Base),
      inherit.aes = FALSE, size = 4, fontface = "bold"
    )
  }

  return(p)
}
