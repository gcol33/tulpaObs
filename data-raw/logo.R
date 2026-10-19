args <- commandArgs(trailingOnly = TRUE)
out <- if (length(args)) args[1] else "dev_notes/hex-icons"
dir.create(out, recursive = TRUE, showWarnings = FALSE)
e <- new.env()
data("crossbill", package = "unmarked", envir = e)
d <- e$crossbill
paper <- "#F5F6F8"
blue <- "#2E6F9E"
teal <- "#2A7A7A"
orange <- "#9E5E0F"
hex <- function(r) {
  a <- seq(-pi / 2, length.out = 6, by = pi / 3)
  paste(sprintf("%.3f,%.3f", 260 + r * cos(a), 300 + r * sin(a)), collapse = " ")
}
circle <- function(x, y, r, fill, stroke = "none", width = 0) {
  sprintf('<circle cx="%.3f" cy="%.3f" r="%.3f" fill="%s" stroke="%s" stroke-width="%.2f"/>',
          x, y, r, fill, stroke, width)
}
path <- function(x, y, color, width = 3, fill = "none") {
  sprintf('<path d="M %s" fill="%s" stroke="%s" stroke-width="%.2f" stroke-linejoin="round" stroke-linecap="round"/>',
          paste(sprintf("%.3f %.3f", x, y), collapse = " L "), fill, color, width)
}
background <- function(kind) {
  if (kind == "radial") {
    paste(vapply(seq(0, 324, 36), function(a) {
      t <- a * pi / 180
      path(c(260, 260 + 600 * cos(t)), c(246, 246 + 600 * sin(t)), "#DDE3E8", 2)
    }, ""), collapse = "")
  } else {
    paste(vapply(seq(-1800, 1800, 36), function(x) {
      if (kind == "diagonal") path(c(x, x + 1800), c(-600, 1200), "#DDE3E8", 2)
      else path(c(-60, 580), c(x, x), "#DDE3E8", 2)
    }, ""), collapse = "")
  }
}
write_logo <- function(id, picture, kind) {
  bg <- background(kind)
  animation <- if (kind == "radial") {
    '<animateTransform attributeName="transform" type="rotate" from="0 260 246" to="36 260 246" dur="12s" repeatCount="indefinite"/>'
  } else {
    '<animateTransform attributeName="transform" type="translate" from="0 0" to="0 36" dur="12s" repeatCount="indefinite"/>'
  }
  svg <- paste0('<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 520 600" role="img" aria-label="tulpaObs logo"><title>tulpaObs logo</title>',
    '<style>.still{display:none}@media(prefers-reduced-motion:reduce){.moving{display:none}.still{display:inline}}</style>',
    '<defs><filter id="dot-shadow" x="-80%" y="-80%" width="260%" height="260%"><feGaussianBlur stdDeviation="4"/></filter><clipPath id="hex"><polygon points="', hex(288), '"/></clipPath>',
    '<linearGradient id="fade" x2="0" y2="1"><stop offset="0.58" stop-color="white"/><stop offset="0.76" stop-color="black"/></linearGradient>',
    '<mask id="soft"><rect width="520" height="600" fill="url(#fade)"/></mask></defs>',
    '<polygon fill="', paper, '" points="', hex(288), '"/>',
    '<g clip-path="url(#hex)"><g mask="url(#soft)"><g class="moving">', animation, bg,
    '</g><g class="still">', bg, '</g></g>', picture, '</g>',
    '<text x="260" y="472" text-anchor="middle" font-family="Georgia, Times New Roman, serif" font-weight="bold" font-style="italic" font-size="76" fill="', blue, '" stroke="', paper, '" stroke-width="10" paint-order="stroke">tulpaObs</text>',
    '<polygon points="', hex(291), '" fill="none" stroke="#2B3A4A" stroke-width="18" stroke-linejoin="round"/></svg>')
  stopifnot(!grepl("NaN|Inf", svg), nchar(svg, type = "bytes") < 200000)
  writeLines(svg, file.path(out, paste0(id, ".svg")))
}

# Consecutive complete survey histories containing at least one detection.
y <- as.matrix(d[, c("det071", "det072", "det073")])
selected <- which(rowSums(y, na.rm = TRUE) > 0 & rowSums(is.na(y)) == 0)[1:7]
write.csv(data.frame(id = d$id[selected], y[selected, ]), file.path(out, "histories.csv"), row.names = FALSE)
history_marks <- function(layout, radius = 11, stroke = 3, track = 2,
                          colors = rep(teal, 3), depth = FALSE, rows = selected) {
  marks <- character()
  u <- seq(0, 1, length.out = 25)
  basis <- cbind((1 - u) * (1 - 2 * u), 4 * u * (1 - u), u * (2 * u - 1))
  for (i in seq_along(rows)) {
    xy <- layout(i)
    if (depth) {
      marks <- c(marks, sprintf('<rect x="%.2f" y="145" width="34" height="210" rx="17" fill="%s" fill-opacity="0.045"/>', xy$x[1] - 17, teal))
      for (j in 1:3) {
        if (y[rows[i], j] == 1) marks <- c(marks,
          sprintf('<g opacity="0.2" filter="url(#dot-shadow)">%s</g>',
                  circle(xy$x[j], xy$y[j] + 5, radius + 2, colors[j])))
      }
    }
    marks <- c(marks, path(as.numeric(basis %*% xy$x), as.numeric(basis %*% xy$y), "#D3DDDF", track))
    for (j in 1:3) marks <- c(marks, circle(xy$x[j], xy$y[j], radius * xy$scale,
      if (y[rows[i], j] == 1) colors[j] else paper, colors[j], stroke * xy$scale))
  }
  paste(marks, collapse = "")
}
flat_history <- function(i) list(x = c(152, 260, 368), y = rep(134 + (i - 1) * 36, 3), scale = 1)
write_logo("01-histories", history_marks(flat_history), "diagonal")

history_variants <- list(
  "01a-bold" = list(kind = "horizontal", radius = 15, stroke = 3.5, track = 3,
    layout = function(i) list(x = c(146, 260, 374), y = rep(126 + (i - 1) * 40, 3), scale = 1)),
  "01b-curved" = list(kind = "diagonal", radius = 14, stroke = 3.2, track = 2.5,
    layout = function(i) {
      t <- (i - 1) / 6
      halfwidth <- 154 * (1 - 0.32 * t^2)
      yy <- 132 + 236 * t
      list(x = 260 + c(-halfwidth, 0, halfwidth), y = yy - c(1, 0, 1) * 28 * t^2,
           scale = 1 - 0.22 * t)
    }),
  "01c-columns" = list(kind = "horizontal", radius = 13, stroke = 3.2, track = 2.5,
    layout = function(i) list(x = rep(104 + (i - 1) * 52, 3), y = c(164, 250, 336), scale = 1))
)
for (id in names(history_variants)) {
  v <- history_variants[[id]]
  write_logo(id, history_marks(v$layout, v$radius, v$stroke, v$track), v$kind)
}
v <- history_variants[["01c-columns"]]
eligible <- which(rowSums(y, na.rm = TRUE) > 0 & rowSums(is.na(y)) == 0)
patterns <- apply(y[eligible, , drop = FALSE], 1, paste0, collapse = "")
balanced <- eligible[!duplicated(patterns)]
stopifnot(length(balanced) == 7, all(colSums(y[balanced, ]) == 4))
write.csv(data.frame(id = d$id[balanced], y[balanced, ]),
          file.path(out, "layered-histories.csv"), row.names = FALSE)
write_logo("01d-layered", history_marks(v$layout, radius = 14, stroke = 3,
  track = 2.5, colors = c(teal, "#7A4F9E", blue), depth = TRUE, rows = balanced), "horizontal")

# Elevation-bin detection frequencies use the 2007 repeated surveys.
bin <- cut(d$ele, breaks = seq(0, 3000, 300), include.lowest = TRUE)
complete <- rowSums(!is.na(y)) > 0
freq <- tapply(rowSums(y, na.rm = TRUE)[complete] > 0, bin[complete], mean)
mid <- seq(150, 2850, 300)
ok <- is.finite(freq)
xx <- 92 + mid[ok] / 3000 * 336
yy <- 354 - as.numeric(freq[ok]) * 220
profile <- c(path(c(xx[1], xx, tail(xx, 1)), c(368, yy, 368), "none", 0, "#DFEAF1"),
  path(xx, yy, blue, 7))
for (i in seq_along(xx)) profile <- c(profile, circle(xx[i], yy[i], 9, paper, blue, 4))
write.csv(data.frame(elevation_midpoint = mid, detection_frequency = as.numeric(freq)),
          file.path(out, "elevation.csv"), row.names = FALSE)
write_logo("02-elevation", paste(profile, collapse = ""), "horizontal")

# Each ring records the annual fraction of surveyed sites with a detection.
years <- c("99", sprintf("%02d", 0:7))
rates <- vapply(years, function(year) {
  z <- as.matrix(d[, paste0("det", year, 1:3)])
  surveyed <- rowSums(!is.na(z)) > 0
  mean(rowSums(z[surveyed, , drop = FALSE], na.rm = TRUE) > 0)
}, 0)
arcs <- character()
for (i in seq_along(rates)) {
  r <- 35 + (i - 1) * 15
  a <- seq(-pi / 2, 3 * pi / 2, length.out = 121)
  arcs <- c(arcs, path(260 + r * cos(a), 246 + r * sin(a), "#DDE3E8", 7))
  a <- seq(-pi / 2, -pi / 2 + 2 * pi * rates[i], length.out = 61)
  arcs <- c(arcs, path(260 + r * cos(a), 246 + r * sin(a), orange, 8))
}
write.csv(data.frame(year = 1999:2007, detected_fraction = rates), file.path(out, "annual.csv"), row.names = FALSE)
write_logo("03-annual", paste(arcs, collapse = ""), "radial")
cat("Built logos from unmarked::crossbill in", out, "\n")
