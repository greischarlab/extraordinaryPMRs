# Code to generate Fig. 4 from Matlab output

# to use Arial fonts in figures:
library("extrafont")
#font_import() # takes awhile only run when needed
#font_import(recursive=FALSE) # takes awhile only run when needed
quartzFonts(avenir=c("Avenir Book", "Avenir Black", "Avenir Book Oblique", "Avenir Black Oblique"),Arial=c("Arial Book", "Arial Black", "Arial Book Oblique", "Arial Black Oblique"))
loadfonts(device = "postscript")

obsPMRShort = read.csv("PMRMaxObs12m.csv", header=F)
regPMRShort = read.csv("PMRFromReg12m.csv", header=F)
obsPMRLong = read.csv("PMRMaxObs6h.csv", header=F)
regPMRLong = read.csv("PMRFromReg6h.csv", header=F)
obsGuitarFlipped = t(read.csv("PMRHeatmapMaxObs.csv", header=F))
obsGuitar = obsGuitarFlipped[,length(obsGuitarFlipped[1,]):1]
regGuitarFlipped = t(read.csv("PMRHeatmapFromReg.csv", header=F))
regGuitar = regGuitarFlipped[,length(regGuitarFlipped[1,]):1]

# color bar that doesn't make a new postscript
color.bar <- function(lut, min, max=-min, nticks=11, ticks=seq(min, max, len=nticks), 
                          labs = ticks, title='', cexVal=1, tickOffset = NA, colWhite = 1) {
  #scale = (length(lut))/(max-min)
  
  plot(c(min,max), c(0,6), type='n', bty='n', xaxt='n', xlab='', yaxt='n', ylab='', 
       main=title, xaxs = "i", yaxs = "i")
  for (i in 1:(length(lut))) {
    x0 = ticks[i]
    x1 = ticks[(i+1)]
    rect(x0,0,x1,2, col=lut[i], border=NA)
    }
  if(is.na(tickOffset)){
    axis(1, at=ticks, labels=labs, las=1, cex.axis=cexVal, tcl = -0.1, mgp = c(3,0.1,0))
  }
  if(!is.na(tickOffset)){
    axis(1, at=ticks, labels=NA, las=1, cex.axis=cexVal, tcl = 0.1, mgp = c(3,0.1,0))
    text(x = ticks[1:length(labs)]+tickOffset, y = rep(1,length(labs)), labs, 
         col = c(rep("white", colWhite), rep("black", length(labs)-colWhite)), cex = cexVal)
    # axis(1, at = ticks[1:length(labs)]+tickOffset, las = 1, labels = labs, tcl = 0, 
    #      cex.axis = cexVal, mgp = c(3,0.1,0))
  }
  box()
}

# modified legend function that allows line end to be altered
legendMod <- function (x, y = NULL, legend, fill = NULL, col = par("col"), 
                       border = "black", lty, lwd, pch, angle = 45, density = NULL, 
                       bty = "o", bg = par("bg"), box.lwd = par("lwd"), box.lty = par("lty"), 
                       box.col = par("fg"), pt.bg = NA, cex = 1, pt.cex = cex, pt.lwd = lwd, 
                       xjust = 0, yjust = 1, x.intersp = 1, y.intersp = 1, adj = c(0, 0.5), 
                       text.width = NULL, text.col = par("col"), text.font = NULL, 
                       merge = do.lines && has.pch, trace = FALSE, plot = TRUE, 
                       ncol = 1, horiz = FALSE, title = NULL, inset = 0, xpd, title.col = text.col[1], 
                       title.adj = 0.5, title.cex = cex[1], title.font = text.font[1], 
                       seg.len = 2, seg.lend = 2) 
{
  if (missing(legend) && !missing(y) && (is.character(y) || 
                                         is.expression(y))) {
    legend <- y
    y <- NULL
  }
  mfill <- !missing(fill) || !missing(density)
  if (!missing(xpd)) {
    op <- par("xpd")
    on.exit(par(xpd = op))
    par(xpd = xpd)
  }
  text.font <- if (is.null(text.font)) 
    par("font")
  else text.font
  title <- as.graphicsAnnot(title)
  if (length(title) > 1) 
    stop("invalid 'title'")
  legend <- as.graphicsAnnot(legend)
  if (any(sapply(legend, is.language))) 
    legend <- as.expression(legend)
  n.leg <- length(legend)
  if (n.leg == 0) 
    stop("'legend' is of length 0")
  auto <- if (is.character(x)) 
    match.arg(x, c("bottomright", "bottom", "bottomleft", 
                   "left", "topleft", "top", "topright", "right", "center"))
  else NA
  if (is.na(auto)) {
    xy <- xy.coords(x, y, setLab = FALSE)
    x <- xy$x
    y <- xy$y
    nx <- length(x)
    if (nx < 1 || nx > 2) 
      stop("invalid coordinate lengths")
  }
  else nx <- 0
  reverse.xaxis <- par("xaxp")[1] > par("xaxp")[2]
  reverse.yaxis <- par("yaxp")[1] > par("yaxp")[2]
  xlog <- par("xlog")
  ylog <- par("ylog")
  cex <- rep(cex, length.out = n.leg)
  x.intersp <- rep(x.intersp, length.out = n.leg)
  seg.len <- rep(seg.len, length.out = n.leg)
  rect2 <- function(left, top, dx, dy, density = NULL, angle, 
                    ...) {
    r <- left + dx
    if (xlog) {
      left <- 10^left
      r <- 10^r
    }
    b <- top - dy
    if (ylog) {
      top <- 10^top
      b <- 10^b
    }
    rect(left, top, r, b, angle = angle, density = density, 
         ...)
  }
  segments2 <- function(x1, y1, dx, dy, ...) {
    x2 <- x1 + dx
    if (xlog) {
      x1 <- 10^x1
      x2 <- 10^x2
    }
    y2 <- y1 + dy
    if (ylog) {
      y1 <- 10^y1
      y2 <- 10^y2
    }
    segments(x1, y1, x2, y2, ..., lend = seg.lend)
  }
  points2 <- function(x, y, ...) {
    if (xlog) 
      x <- 10^x
    if (ylog) 
      y <- 10^y
    points(x, y, ...)
  }
  text2 <- function(x, y, ...) {
    if (xlog) 
      x <- 10^x
    if (ylog) 
      y <- 10^y
    text(x, y, ...)
  }
  colwise <- function(x, n, ncol, n.legpercol, fun, reverse = FALSE) {
    xmat <- matrix(c(rep(x, length.out = n), rep(0L, n.legpercol * 
                                                   ncol - n)), ncol = ncol)
    res <- apply(xmat, 2, fun)
    res[res == 0L] <- max(res)
    if (reverse) 
      -res
    else res
  }
  rowwise <- function(x, n, ncol, n.legpercol, fun, reverse = FALSE) {
    xmat <- matrix(c(rep(x, length.out = n), rep(0L, n.legpercol * 
                                                   ncol - n)), ncol = ncol)
    res <- apply(xmat, 1, fun)
    if (reverse) 
      -res
    else res
  }
  if (trace) {
    catn <- function(...) do.call(cat, c(lapply(list(...), 
                                                formatC), "\n"))
    fv <- function(...) paste(vapply(lapply(list(...), formatC), 
                                     paste, collapse = ",", ""), collapse = ", ")
  }
  n.legpercol <- if (horiz) {
    if (ncol != 1) 
      warning(gettextf("horizontal specification overrides: Number of columns := %d", 
                       n.leg), domain = NA)
    ncol <- n.leg
    1
  }
  else ceiling(n.leg/ncol)
  Cex <- cex * par("cex")
  if (is.null(text.width)) 
    text.width <- max(abs(mapply(strwidth, legend, cex = cex, 
                                 font = text.font, MoreArgs = list(units = "user"))))
  else if ((length(text.width) > 1L && any(is.na(text.width))) || 
           (all(!is.na(text.width)) && (!is.numeric(text.width) || 
                                        any(text.width < 0)))) 
    stop("'text.width' must be numeric, >= 0, or a scalar NA")
  if (auto.text.width <- all(is.na(text.width))) {
    text.width <- abs(mapply(strwidth, legend, cex = cex, 
                             font = text.font, MoreArgs = list(units = "user")))
    ncol <- ceiling(n.leg/n.legpercol)
  }
  xyc <- xyinch(par("cin"), warn.log = FALSE)
  xc <- Cex * xyc[1L]
  yc <- Cex * xyc[2L]
  if (any(xc < 0)) 
    text.width <- -text.width
  xchar <- xc
  xextra <- 0
  y.intersp <- rep(y.intersp, length.out = n.legpercol)
  yextra <- rowwise(yc, n = n.leg, ncol = ncol, n.legpercol = n.legpercol, 
                    fun = function(x) max(abs(x)), reverse = reverse.yaxis) * 
    (y.intersp - 1)
  ymax <- sign(yc[1]) * max(abs(yc)) * max(1, mapply(strheight, 
                                                     legend, cex = cex, font = text.font, MoreArgs = list(units = "user"))/yc)
  ychar <- yextra + ymax
  ymaxtitle <- title.cex * par("cex") * xyc[2L] * max(1, strheight(title, 
                                                                   cex = title.cex, font = title.font, units = "user")/(title.cex * 
                                                                                                                          par("cex") * xyc[2L]))
  ychartitle <- yextra[1] + ymaxtitle
  if (trace) 
    catn("  xchar=", fv(xchar), "; (yextra, ychar)=", fv(yextra, 
                                                         ychar))
  if (mfill) {
    xbox <- xc * 0.8
    ybox <- yc * 0.5
    dx.fill <- max(xbox)
  }
  do.lines <- (!missing(lty) && (is.character(lty) || any(lty > 
                                                            0))) || !missing(lwd)
  has.pch <- !missing(pch) && length(pch) > 0
  if (do.lines) {
    x.off <- if (merge) 
      -0.7
    else 0
  }
  else if (merge) 
    warning("'merge = TRUE' has no effect when no line segments are drawn")
  if (has.pch) {
    if (is.character(pch) && !is.na(pch[1L]) && nchar(pch[1L], 
                                                      type = "c") > 1) {
      if (length(pch) > 1) 
        warning("not using pch[2..] since pch[1L] has multiple chars")
      np <- nchar(pch[1L], type = "c")
      pch <- substr(rep.int(pch[1L], np), 1L:np, 1L:np)
    }
    if (!is.character(pch)) 
      pch <- as.integer(pch)
  }
  if (is.na(auto)) {
    if (xlog) 
      x <- log10(x)
    if (ylog) 
      y <- log10(y)
  }
  if (nx == 2) {
    x <- sort(x)
    y <- sort(y)
    left <- x[1L]
    top <- y[2L]
    w <- diff(x)
    h <- diff(y)
    w0 <- w/ncol
    x <- mean(x)
    y <- mean(y)
    if (missing(xjust)) 
      xjust <- 0.5
    if (missing(yjust)) 
      yjust <- 0.5
  }
  else {
    yc <- rowwise(yc, n.leg, ncol, n.legpercol, fun = function(x) max(abs(x)), 
                  reverse = reverse.yaxis)
    h <- sum(ychar) + yc[length(yc)] + (!is.null(title)) * 
      ychartitle
    xch1 <- colwise(xchar, n.leg, ncol, n.legpercol, fun = function(x) max(abs(x)), 
                    reverse = reverse.xaxis)
    x.interspCol <- colwise(x.intersp, n.leg, ncol, n.legpercol, 
                            fun = max)
    seg.lenCol <- colwise(seg.len, n.leg, ncol, n.legpercol, 
                          fun = max)
    text.width <- colwise(text.width, n = if (auto.text.width) 
      n.leg
      else ncol, ncol, n.legpercol = if (auto.text.width) 
        n.legpercol
      else 1, fun = function(x) max(abs(x)), reverse = reverse.xaxis)
    w0 <- text.width + (x.interspCol + 1) * xch1
    if (mfill) 
      w0 <- w0 + dx.fill
    if (do.lines) 
      w0 <- w0 + (seg.lenCol + x.off) * xch1
    w <- sum(w0) + 0.5 * xch1[ncol]
    if (!is.null(title) && (abs(tw <- strwidth(title, units = "user", 
                                               cex = title.cex, font = title.font) + 0.5 * title.cex * 
                                par("cex") * xyc[1L])) > abs(w)) {
      xextra <- (tw - w)/2
      w <- tw
    }
    if (is.na(auto)) {
      left <- x - xjust * w
      top <- y + (1 - yjust) * h
    }
    else {
      usr <- par("usr")
      inset <- rep_len(inset, 2)
      insetx <- inset[1L] * (usr[2L] - usr[1L])
      left <- switch(auto, bottomright = , topright = , 
                     right = usr[2L] - w - insetx, bottomleft = , 
                     left = , topleft = usr[1L] + insetx, bottom = , 
                     top = , center = (usr[1L] + usr[2L] - w)/2)
      insety <- inset[2L] * (usr[4L] - usr[3L])
      top <- switch(auto, bottomright = , bottom = , bottomleft = usr[3L] + 
                      h + insety, topleft = , top = , topright = usr[4L] - 
                      insety, left = , right = , center = (usr[3L] + 
                                                             usr[4L] + h)/2)
    }
  }
  if (plot && bty != "n") {
    if (trace) 
      catn("  rect2(", left, ",", top, ", w=", w, ", h=", 
           h, ", ...)", sep = "")
    rect2(left, top, dx = w, dy = h, col = bg, density = NULL, 
          lwd = box.lwd, lty = box.lty, border = box.col)
  }
  xt <- left + xc + xextra + rep(c(0, cumsum(w0))[1L:ncol], 
                                 each = n.legpercol, length.out = n.leg)
  topspace <- 0.5 * ymax + (!is.null(title)) * ychartitle
  yt <- top - topspace - cumsum((c(0, ychar)/2 + c(ychar, 0)/2)[1L:n.legpercol])
  yt <- rep(yt, length.out = n.leg)
  if (mfill) {
    if (plot) {
      if (!is.null(fill)) 
        fill <- rep_len(fill, n.leg)
      rect2(left = xt, top = yt + ybox/2, dx = xbox, dy = ybox, 
            col = fill, density = density, angle = angle, 
            border = border)
    }
    xt <- xt + dx.fill
  }
  if (plot && (has.pch || do.lines)) 
    col <- rep_len(col, n.leg)
  if (missing(lwd) || is.null(lwd)) 
    lwd <- par("lwd")
  if (do.lines) {
    if (missing(lty) || is.null(lty)) 
      lty <- 1
    lty <- rep_len(lty, n.leg)
    lwd <- rep_len(lwd, n.leg)
    ok.l <- !is.na(lty) & (is.character(lty) | lty > 0) & 
      !is.na(lwd)
    if (trace) 
      catn("  segments2(", xt[ok.l] + x.off * xchar[ok.l], 
           ",", yt[ok.l], ", dx=", (seg.len * xchar)[ok.l], 
           ", dy=0, ...)")
    if (plot) 
      segments2(xt[ok.l] + x.off * xchar[ok.l], yt[ok.l], 
                dx = (seg.len * xchar)[ok.l], dy = 0, lty = lty[ok.l], 
                lwd = lwd[ok.l], col = col[ok.l])
    xt <- xt + (seg.len + x.off) * xchar
  }
  if (has.pch) {
    pch <- rep_len(pch, n.leg)
    pt.bg <- rep_len(pt.bg, n.leg)
    pt.cex <- rep_len(pt.cex, n.leg)
    pt.lwd <- rep_len(pt.lwd, n.leg)
    ok <- !is.na(pch)
    if (!is.character(pch)) {
      ok <- ok & (pch >= 0 | pch <= -32)
    }
    else {
      ok <- ok & nzchar(pch)
    }
    x1 <- (if (merge && do.lines) 
      xt - (seg.len/2) * xchar
      else xt)[ok]
    y1 <- yt[ok]
    if (trace) 
      catn("  points2(", x1, ",", y1, ", pch=", pch[ok], 
           ", ...)")
    if (plot) 
      points2(x1, y1, pch = pch[ok], col = col[ok], cex = pt.cex[ok], 
              bg = pt.bg[ok], lwd = pt.lwd[ok])
  }
  xt <- xt + x.intersp * xc
  if (plot) {
    if (!is.null(title)) 
      text2(left + w * title.adj, top - ymaxtitle, labels = title, 
            adj = c(title.adj, 0), cex = title.cex, col = title.col, 
            font = title.font)
    text2(xt, yt, labels = legend, adj = adj, cex = cex, 
          col = text.col, font = text.font)
  }
  invisible(list(rect = list(w = w, h = h, left = left, top = top), 
                 text = list(x = xt, y = yt)))
}


medAge = c(0,6,12,18,24,30,36,42) # x values throughout
truePMR = seq(2,32,by=2) # y values for heatmap
xOffset = 0.4
lwdVal = 4
midline = 0.5
xAxisLine = 0.8
yAxisLine = 1.2
xPlotLabLine = -0.6
yPlotLabLine = 2
defaultMarLeft = c(2,3,0.5,0)
defaultMarRight = c(2,2,0.5,1)
keyMarLeft = c(0,4,0,1)
keyMarRight = c(0,3,0,2)

obsColorPalette = colorRampPalette(c("darkred","white"))
regColorPalette = colorRampPalette(c("cornflowerblue","white"))

cexText = 0.8
titleExp = 1
lowPMRObsCol = obsColorPalette(3)[1]
hiPMRObsCol = obsColorPalette(3)[2]
lowPMRRegCol = regColorPalette(3)[1]
hiPMRRegCol = regColorPalette(3)[2]

heatMapXFrac = 1
keyXFrac = 0
keyYFrac = 0.8
keyUpperYFrac = 0.92

postscript("Fig4.eps", width = 6.5, height = 4.25, paper="special", horizontal=F, family="ArialMT")

par(fig = c(0,midline,0.5,1), mar = defaultMarLeft, bty='n', xpd = NA)

plot(range(medAge), range(obsPMRShort), type='n', xlab = "", ylab = "",
     xaxt='n', yaxt = 'n', ylim = log10(c(5,1000)))
axis(1, at = seq(0, 42, by = 6), tcl = -0.1, mgp = c(3,0.05,0), cex.axis = cexText)
mtext("Initial median parasite age", side = 1, line = xAxisLine, cex = cexText)
axis(2, tcl = -0.1, at = log10(c(10,100,1000)), 
     labels = c(10,expression(10^2),expression(10^3)),
     mgp = c(3,0.2,0), las=2, cex.axis = cexText)
mtext("Estimated PMR", side = 2, line = yAxisLine, cex = cexText)
mtext(expression(bold("12 min. sampling window")), side = 2, line = yPlotLabLine, cex = titleExp*cexText)
mtext(expression(phantom("g")~bold("Maximum observed")~phantom("g")), side=3, 
      line = xPlotLabLine, cex = titleExp*cexText)
segments(x0 = medAge-xOffset, y0 = log10(obsPMRShort[,1]), y1 = log10(obsPMRShort[,2]), 
         col = lowPMRObsCol, lwd = lwdVal, lend = 2)
segments(x0 = medAge+xOffset, y0 = log10(obsPMRShort[,3]), y1 = log10(obsPMRShort[,4]), col = hiPMRObsCol, 
         lwd = lwdVal, lend = 2)
text(x=0,y = 3, "(A)", adj=0, cex = cexText)
legendMod(x = 22, y = 2.8, legend = c(6,12), title = "True PMR", col = c(lowPMRObsCol, hiPMRObsCol), 
          lwd = lwdVal, bty='n', cex = cexText)

par(new=T,fig = c(midline,1,0.5,1), mar = defaultMarRight, bty='n')
plot(range(medAge), range(regPMRShort), type='n', xlab = "", ylab = "",
     xaxt='n', yaxt = 'n', ylim = log10(c(5,1000)))
axis(1, at = seq(0, 42, by = 6), tcl = -0.1, mgp = c(3,0.05,0), cex.axis = cexText)
axis(2, tcl = -0.1, at = log10(c(10,100,1000)), 
     labels = c(10,expression(10^2),expression(10^3)),
     mgp = c(3,0.2,0), las=2, cex.axis = cexText)
mtext("Initial median parasite age", side = 1, line = xAxisLine, cex = cexText)
mtext(expression(bold("Regression")), side=3, line = xPlotLabLine, cex = titleExp*cexText)
segments(x0 = medAge-xOffset, y0 = log10(regPMRShort[,1]), y1 = log10(regPMRShort[,2]), 
         col = lowPMRRegCol, lwd = lwdVal, lend = 2)
segments(x0 = medAge+xOffset, y0 = log10(regPMRShort[,3]), y1 = log10(regPMRShort[,4]), 
         col = hiPMRRegCol, lwd = lwdVal, lend = 2)
text(x=0,y=3,"(B)",adj=0, cex = cexText)
legendMod(x = 22, y = 2.8, legend = c(6,12), title = "True PMR", col = c(lowPMRRegCol, hiPMRRegCol), 
          lwd = lwdVal, bty='n', cex = cexText)

par(new=T,fig = c(0,midline, 0, 0.5), mar = defaultMarLeft, bty='n')
plot(range(medAge), range(obsPMRLong), type='n', xlab = "", ylab = "",
     xaxt='n', yaxt = 'n', ylim = log10(c(5,1000)))
axis(1, at = seq(0, 42, by = 6), tcl = -0.1, mgp = c(3,0.05,0), cex.axis = cexText)
mtext("Initial median parasite age", side = 1, line = xAxisLine, cex = cexText)
axis(2, tcl = -0.1, at = log10(c(10,100,1000)), 
     labels = c(10,expression(10^2),expression(10^3)),
     mgp = c(3,0.2,0), las=2, cex.axis = cexText)
mtext("Estimated PMR", side = 2, line = yAxisLine, cex = cexText)
mtext(expression(bold("6 hour sampling window")), side = 2, line = yPlotLabLine, cex = titleExp*cexText)
segments(x0 = medAge-xOffset, y0 = log10(obsPMRLong[,1]), y1 = log10(obsPMRLong[,2]), 
         col = lowPMRObsCol, lwd = lwdVal, lend = 2)
segments(x0 = medAge+xOffset, y0 = log10(obsPMRLong[,3]), y1 = log10(obsPMRLong[,4]), 
         col = hiPMRObsCol, lwd = lwdVal, lend = 2)
text(x=0,y = 3, "(C)", adj=0, cex = cexText)

par(new=T,fig = c(midline,1, 0, 0.5), mar = defaultMarRight, bty='n')
plot(range(medAge), range(regPMRLong), type='n', xlab = "", ylab = "",
     xaxt='n', yaxt = 'n', ylim = log10(c(5,1000)))
axis(1, at = seq(0, 42, by = 6), tcl = -0.1, mgp = c(3,0.05,0),cex.axis = cexText)
axis(2, tcl = -0.1, at = log10(c(10,100,1000)), 
     labels = c(10,expression(10^2),expression(10^3)),
     mgp = c(3,0.2,0), las=2, cex.axis = cexText)
mtext("Initial median parasite age", side = 1, line = xAxisLine, cex = cexText)
segments(x0 = medAge-xOffset, y0 = log10(regPMRLong[,1]), y1 = log10(regPMRLong[,2]), 
         col = lowPMRRegCol, lwd = lwdVal, lend = 2)
segments(x0 = medAge+xOffset, y0 = log10(regPMRLong[,3]), y1 = log10(regPMRLong[,4]), 
         col = hiPMRRegCol, lwd = lwdVal, lend = 2)
text(x=0,y=3,"(D)",adj=0, cex = cexText)

dev.off()

postscript("FigS8.eps", width = 6.5, height = 3.25, paper="special", horizontal=F, family="ArialMT")

par(fig = c(0, heatMapXFrac*midline, 0, 1), 
    mar = defaultMarLeft, bty='n', xpd = NA)

# need one fewer color than break
breakValues = seq(0,6.2,by=0.2) - 0.1
# make one too many colors so we can drop the last one (white)
colValues = c("black", obsColorPalette(length(breakValues)-1))
image(x = medAge, y = truePMR, z = obsGuitar, xlab = "", ylab = "",
      breaks = breakValues, col = colValues[1:(length(colValues)-1)], 
      xaxt = 'n', yaxt = 'n', ylim = c(1,43))
axis(1, at = seq(0, 42, by = 6), tcl = -0.1, mgp = c(3,0.05,0), cex.axis = cexText)
mtext("Initial median parasite age", side = 1, line = xAxisLine, cex = cexText)
axis(2, tcl = -0.1, at = truePMR, labels = c(2, "", 6, "", 10, "", 14, "", 18, "", 22, "", 26, "", 30, ""),
     mgp = c(3,0.2,0), las=2, cex.axis = cexText)
mtext("True PMR compared to PMR = 6", side = 2, line = yAxisLine, cex = cexText)
mtext(expression(bold("Sampling window required")), side = 2, 
      line = yPlotLabLine, cex = titleExp*cexText)
text(x = 0, y = 42, "(A)", cex = cexText)
for (val in 1:length(truePMR)){
  colVals = rep("black", length(medAge))
  textVals = signif(obsGuitar[,val], digits = 2)
  # three special cases:
  if(any(textVals<=2, na.rm = T)){colVals[textVals<=2] = "white"}
  if(any(textVals=="0")){
    colVals[textVals=="0"] = "white"
    textVals[textVals=="0"] = ""
    }
  if(any(is.nan(textVals))){textVals[is.nan(textVals)] <- ""}
  text(x=medAge,y = rep(truePMR[val],length(medAge)),textVals,cex=0.7,col = colVals)
}
text(x = 21, y = 6, "NA, identical PMRs", col = "white", cex = 0.7)
text(x = 22, y = 32.8, "Largest window\nwith no overlap\nexceeds 6 hrs", adj = c(0,1), cex = 0.7)

# vector of values to plot on color bar (excepting the NaN values)
toPlot = which(signif(breakValues+0.1, digits=2) %in% signif(c(0:6), digits = 2))

breakValLabs = c(breakValues[toPlot]+0.1, 7:8)
breakValLabs[1] = "NA"
breakValLabs[(length(breakValLabs)-1)] = ">6"
breakValLabs[length(breakValLabs)] = "(hrs)"

par(new=T,fig = c(keyXFrac*midline,midline, 
                  keyYFrac, keyUpperYFrac), 
    mar = keyMarLeft)
color.bar(c(colValues[toPlot],"white"), min = 0, max=8, 
          ticks=0:8, labs = breakValLabs, tickOffset = 0.5, cex = 0.7*cexText, colWhite = 3)
mtext("Largest sampling window with no overlap", line = -1, cex = cexText)

par(new=T,fig = c(midline,midline + heatMapXFrac*midline, 0, 1), 
    mar = defaultMarRight, bty='n')

# need one fewer color than break
breakValues = seq(0,6.2,by=0.2) - 0.1
# make one too many colors so we can drop the last one (white)
colValues = c("black", regColorPalette(length(breakValues)-1))
image(x = medAge, y = truePMR, z = regGuitar, xlab = "", ylab = "",
      breaks = breakValues, col = colValues[1:(length(colValues)-1)], 
      xaxt = 'n', yaxt = 'n', ylim = c(1,43))
axis(1, at = seq(0, 42, by = 6), tcl = -0.1, mgp = c(3,0.05,0),cex.axis = cexText)
mtext("Initial median parasite age", side = 1, line = xAxisLine, cex = cexText)
axis(2, tcl = -0.1, at = truePMR, labels = c(2, "", 6, "", 10, "", 14, "", 18, "", 22, "", 26, "", 30, ""),
     mgp = c(3,0.2,0), las=2, cex.axis = cexText)
text(x = 0, y = 42, "(B)", cex = cexText)
regGuitarLab = signif(regGuitar, digits = 2);regGuitarLab[is.nan(regGuitarLab)]<-""
for (val in 1:length(truePMR)){
  text(x=medAge,y = rep(truePMR[val],length(medAge)),regGuitarLab[,val],cex=0.7)
}
if(any(regGuitarLab=="0")){
  lessThan = which(regGuitarLab=="0", arr.ind = T)
  for (i in 1:length(lessThan[,1])){
    text(x = medAge[lessThan[i,1]], y = truePMR[lessThan[i,2]], "", col = "white", cex=0.7, adj=c(0.5,0.5))
  }
  }
text(x = 21, y = 6, "NA, identical PMRs", col = "white", cex = 0.7)
text(x = 22, y = 32.8, "Largest window\nwith no overlap\nexceeds 6 hrs", adj = c(0,1), cex = 0.7)

# vector of values to plot on color bar (excepting the NaN values)
toPlot = which(signif(breakValues+0.1, digits=2) %in% signif(c(0:6), digits = 2))

breakValLabs = c(breakValues[toPlot]+0.1, 7:8)
breakValLabs[1] = "NA"
breakValLabs[(length(breakValLabs)-1)] = ">6"
breakValLabs[length(breakValLabs)] = "(hrs)"

par(new=T,fig = c(midline + keyXFrac*midline, 1, keyYFrac, 
                  keyUpperYFrac), mar = keyMarRight)
color.bar(c(colValues[toPlot],"white"), min = 0, max=8, 
          ticks=0:8, labs = breakValLabs, tickOffset = 0.5,cexVal = 0.7*cexText)
mtext("Largest sampling window with no overlap", line = -1, cex = cexText)

dev.off()
