####################################################################
# This script is for the palaeo-proxies plot                       #                       
# Author: Mathias Bellat                                           #
# Affiliation : Tubingen University                                #
# Creation date : 18/08/2026                                       #
# E-mail: mathias.archaeology@gmail.com                            #
####################################################################

raw <- read.csv("./analysis/data/raw_data/zscores_raw.csv", stringsAsFactors = FALSE)
resampled <- read.csv("./analysis/data/raw_data/zscores_resampled_25yr.csv", stringsAsFactors = FALSE)

series_cols <- c(
  "Jeita_cave", "Soreq_cave", "Shalaii_cave", "KunaBa_cave", "LoNAP514",
  "Lake_Van", "Lake_Zeribar", "Lake_Mirabad",
  "Lake_Urmia_Quercus", "Lake_Van_Quercus", "Lake_Zeribar_Quercus",
  "Lake_Mirabad_Quercus", "Hashilan_Quercus"
)

n_series    <- length(series_cols)
offset_step <- 3      # vertical spacing between trace bands
band_height <- 2.4    # vertical space used by the trace within its band

is_pollen <- function(name) grepl("_Quercus$", name)

palette13 <- c(
  "#1f77b4", "#aec7e8", "#ff7f0e", "#2ca02c", "#d62728",
  "#9467bd", "#8c564b", "#e377c2", "#7f7f7f", "#bcbd22",
  "#17becf", "#9edae5", "#393b79"
)

bp_to_bc <- function(bp) bp - 1950   # astronomical convention (present = 1950 CE)

bp_breaks <- seq(2000, 8000, 1000)
x_max <- 8300   # left edge of the plot (oldest)
x_min <- 1700   # right edge of the plot (most recent), with room for labels

# ---------------------------------------------------------------------------
# Pre-compute, for each series: rescaled trace (native value -> vertical
# band position) and native-unit tick positions
# ---------------------------------------------------------------------------

traces <- vector("list", n_series)
names(traces) <- series_cols

for (i in seq_along(series_cols)) {
  nm  <- series_cols[i]
  sub <- raw[raw$series == nm, ]
  sub <- sub[order(sub$BP), ]
  
  vmin <- min(sub$value, na.rm = TRUE)
  vmax <- max(sub$value, na.rm = TRUE)
  offset <- (i - 1) * offset_step
  
  scale_y <- function(v, vmin_ = vmin, vmax_ = vmax, offset_ = offset) {
    offset_ + (v - vmin_) / (vmax_ - vmin_) * band_height
  }
  
  bp_vec <- sub$BP
  y_vec  <- scale_y(sub$value)
  
  gaps <- diff(bp_vec)
  med_gap <- median(gaps)
  hiatus_threshold <- max(500, 8 * med_gap)
  hiatus_at <- which(gaps > hiatus_threshold)
  
  if (length(hiatus_at) > 0) {
    bp_out <- numeric(0); y_out <- numeric(0)
    prev <- 1
    for (h in hiatus_at) {
      bp_out <- c(bp_out, bp_vec[prev:h], mean(bp_vec[h:(h + 1)]))
      y_out  <- c(y_out,  y_vec[prev:h],  NA)
      prev <- h + 1
    }
    bp_out <- c(bp_out, bp_vec[prev:length(bp_vec)])
    y_out  <- c(y_out,  y_vec[prev:length(bp_vec)])
    bp_vec <- bp_out; y_vec <- y_out
  }
  
  ticks_native <- pretty(c(vmin, vmax), n = 3)
  ticks_native <- ticks_native[ticks_native >= vmin & ticks_native <= vmax]
  
  traces[[nm]] <- list(
    name = nm, bp = bp_vec, y = y_vec,
    offset = offset, vmin = vmin, vmax = vmax,
    ticks_native = ticks_native, ticks_y = scale_y(ticks_native),
    n = nrow(sub), colour = palette13[i],
    hiatuses_yr = if (length(hiatus_at) > 0) gaps[hiatus_at] else numeric(0)
  )
}

# Report any detected hiatuses to the console for the record
for (tr in traces) {
  if (length(tr$hiatuses_yr) > 0) {
    message(sprintf("Hiatus detected in %s: %s year gap(s)",
                    tr$name, paste(round(tr$hiatuses_yr), collapse = ", ")))
  }
}

# ---------------------------------------------------------------------------
# Regional composite (delta-18O proxies only, pollen excluded)
# ---------------------------------------------------------------------------

non_pollen_cols <- series_cols[!is_pollen(series_cols)]
comp_mean   <- rowMeans(resampled[, non_pollen_cols], na.rm = TRUE)
comp_median <- apply(resampled[, non_pollen_cols], 1, median, na.rm = TRUE)
bp_grid     <- resampled$BP

# ---------------------------------------------------------------------------
# Plot
# ---------------------------------------------------------------------------

pdf("./analysis/data/derivated_data/Figure_5_raw.pdf", width = 12, height = 10)

## ---- Panel 1: individual traces, each on its own native-unit axis --------
par(mar = c(3, 2, 4, 6), xaxs = "i")

plot(NA, xlim = c(x_max, x_min), ylim = c(-1, n_series * offset_step + 1),
     xlab = "", ylab = "", axes = FALSE)

for (tr in traces) lines(tr$bp, tr$y, col = tr$colour, lwd = 1)

# Years BC along the top, Years BP along the bottom, of this panel
axis(3, at = bp_breaks, labels = bp_to_bc(bp_breaks), lwd = 1, cex.axis = 0.9)
mtext("Years BC", side = 3, line = 2.4, cex = 0.9)
axis(1, at = bp_breaks, labels = bp_breaks, lwd = 1, cex.axis = 0.9)
mtext("Years BP", side = 1, line = 2.0, cex = 0.9)
title(main = "Individual records at native temporal resolution and native units",
      line = 3.6)
box()

# Trace name + proxy type + n, printed at the end of each line
for (tr in traces) {
  proxy_txt <- if (is_pollen(tr$name)) "% pollen" else "\u03b418O"
  text(x = min(tr$bp) - 60, y = tr$offset + band_height / 2,
       labels = paste0(tr$name, " (", proxy_txt, ", n=", tr$n, ")"),
       col = tr$colour, cex = 0.62, adj = 0)
}

# Small native-unit axis + label for each trace, on the right-hand margin.
# NOTE: axis() draws a full-height vertical line every time it is called, so
# calling it 13 times with a coloured `col` would leave only the LAST
# colour visible on that line (earlier ones painted over). We therefore draw
# the line itself only once (in black, on the first call) and use
# `col.ticks` (R >= 3.6) to colour only the tick marks / labels per trace.
for (i in seq_along(traces)) {
  tr <- traces[[i]]
  axis(4, at = tr$ticks_y, labels = tr$ticks_native, lwd = 1, line = 0.2,
       cex.axis = 0.55, mgp = c(0, 0.4, 0),
       col = if (i == 1) "black" else NA,
       col.ticks = tr$colour, col.axis = tr$colour, las = 1)
  mid_y <- mean(tr$ticks_y)
  if (is_pollen(tr$name)) {
    mtext(side = 4, at = mid_y, text = "% pollen", line = 2.1,
          cex = 0.55, col = tr$colour, las = 1)
  } else {
    mtext(side = 4, at = mid_y, text = expression(delta^{18} * O), line = 2.1,
          cex = 0.55, col = tr$colour, las = 1)
  }
}

## ---- Panel 2: composite stack (delta-18O proxies only) -------------------
par(mar = c(4.5, 4, 3, 6), xaxs = "i")

plot(NA, xlim = c(x_max, x_min),
     ylim = range(c(comp_mean, comp_median), na.rm = TRUE) * 1.1,
     xlab = "Years BP", ylab = "Composite z-score", axes = FALSE)

# Fill the area between the curve and zero, segment by segment, so that
# gaps (NA, where no underlying series has data) never get bridged into a
# spurious polygon -- this is what caused the earlier colour/shape glitches.
draw_area <- function(x, y, col_pos, col_neg) {
  valid <- !is.na(y)
  run_id <- cumsum(c(TRUE, diff(valid) != 0))
  for (id in unique(run_id[valid])) {
    idx <- which(run_id == id & valid)
    if (length(idx) < 2) next
    xi <- x[idx]; yi <- y[idx]
    pos <- pmax(yi, 0); neg <- pmin(yi, 0)
    polygon(c(xi, rev(xi)), c(pos, rep(0, length(pos))),
            col = col_pos, border = NA)
    polygon(c(xi, rev(xi)), c(neg, rep(0, length(neg))),
            col = col_neg, border = NA)
  }
}

draw_area(bp_grid, comp_mean,
          col_pos = adjustcolor("steelblue", alpha.f = 0.25),
          col_neg = adjustcolor("firebrick", alpha.f = 0.25))

lines(bp_grid, comp_mean,   col = "steelblue",  lwd = 1.6)
lines(bp_grid, comp_median, col = "darkorange", lwd = 1.3, lty = 2)
abline(h = 0, col = "black", lwd = 1)

# Years BC along the top, Years BP along the bottom, of this panel too
axis(3, at = bp_breaks, labels = bp_to_bc(bp_breaks), cex.axis = 0.9)
mtext("Years BC", side = 3, line = 2.0, cex = 0.9)
axis(1, at = bp_breaks, cex.axis = 0.9)
axis(2, cex.axis = 0.9, las = 1)
box()

legend("topright", legend = c("Composite (mean)", "Composite (median)"),
       col = c("steelblue", "darkorange"), lty = c(1, 2), lwd = c(1.6, 1.3),
       bty = "n", cex = 0.8)

title(main = "Regional composite stack \u2014 \u03b418O proxies only, pollen excluded",
      line = -1.3)

mtext(paste0(
  "Note: top panel shows raw values at each record's native resolution and units ",
  "(13 series, including pollen).\n",
  "Bottom panel shows the composite z-score from the 8 \u03b418O series only (",
  paste(non_pollen_cols, collapse = ", "), "); pollen proxies excluded."
), side = 1, line = 3.5, cex = 0.55, adj = 0, font = 3)

dev.off()
