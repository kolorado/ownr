# Various functions
# Change log: 
# 2015-07-15: lineplot now accepts factors as the first effect. 

# Stat functions---------
#scaled mass index as a measure for body condition according to Peig 2009
smi  = function (M, L) 
{
  plot(log(M)~log(L))
{
    if (require(smatr)) {
      SMA = sma(log(M)~log(L))
      bSMA = coef(SMA)[2]} 
    else {
      OLS = lm(log(M)~log(L))
      bOLS = coef(OLS)[2]
      r = cor.test(~log(M)+log(L), method = "pearson")$estimate
      #outliers = which(abs(rstandard(ols))>3)
      bSMA = bOLS/r }
  }
  L0 = median(L, na.rm = T)
  SMi = M*((L0/L)^bSMA)
  return(SMi)
}    

# standard error of mean
se <- function(x) {
  sqrt(var(x,na.rm=TRUE)/length(na.omit(x)))
}

#z scores, aka standardized values with 0 mean and unit SD
z = function(x) {(x - mean(x, na.rm = T)) / sd(x, na.rm = T)}


#Repeatability calculation---------- 
# following Lessells & Boag 1994

rpt <- function(x, ...) UseMethod("rpt")


rpt.formula <- function(formula, data)
{
  if (missing(formula) || (length(formula) != 3L) || (length(attr(terms(formula[-2L]),
        "term.labels")) != 1L))
        stop("'formula' missing or incorrect")
  mf <- model.frame(formula = formula, data = data)
  resp.index <- attr(attr(mf, "terms"), "response")
  eff <- factor(mf[[-resp.index]])
  resp <-  mf[[resp.index]]
  #resp <- model.response(mf)   #the same
  model = lm(resp~eff)
  MSa <- anova(model)["Mean Sq"]["eff",]
  MSw <- anova(model)["Mean Sq"]["Residuals",]
  n <- table(eff)
  K <- length(levels(eff))
  M <- sum(n)
  n0 <-  (1/(K-1))* (M - (sum(n^2)/M))
  s2 <- MSw
  s2A <- (MSa-MSw)/n0
  r <- s2A/(s2 + s2A)
  p <- round(anova(model)["Pr(>F)"]["eff",],3)
  F <- round(anova(model)["F value"]["eff",], 3)
  Fdf1 <- anova(model)["Df"]["eff",]
  Fdf2 <- anova(model)["Df"]["Residuals",]
  {if (p == 0) cat("r = ", round(r,3),  ", F", Fdf1, ",", Fdf2, "= ", F, ", p < 0.001", "\n", sep = "")
     else      cat("r = ", round(r,3),  ", F", Fdf1, ",", Fdf2, "= ", F, ", p = ", p, "\n", sep = "")
  }
}


rpt.default <- function (data, response, effect)
{
  m <- match.call(expand.dots = FALSE)
  eff = eval(m$effect, data)
  resp = eval(m$response, data)
  if (is.factor(eff) != TRUE) eff = factor(eff)
  model = lm(resp~eff, data = data)
  MSa <- anova(model)["Mean Sq"]["eff",]
  MSw <- anova(model)["Mean Sq"]["Residuals",]
  n <- table(eff)
  K <- length(levels(eff))
  M <- sum(n)
  n0 <-  (1/(K-1))* (M - (sum(n^2)/M))
  s2 <- MSw
  s2A <- (MSa-MSw)/n0
  r <- s2A/(s2 + s2A)
  p <- round(anova(model)["Pr(>F)"]["eff",],3)
  F <- round(anova(model)["F value"]["eff",], 3)
  Fdf1 <- anova(model)["Df"]["eff",]
  Fdf2 <- anova(model)["Df"]["Residuals",]
  {if (p == 0) cat("r = ", round(r,3),  ", F", Fdf1, ",", Fdf2, "= ", F, ", p < 0.001", "\n", sep = "")
     else      cat("r = ", round(r,3),  ", F", Fdf1, ",", Fdf2, "= ", F, ", p = ", p, "\n", sep = "")
  }
}


#--------------Errorbar function to display mean+CI----------------------------
errorbars<-function (response, factor1, factor2, error.bars = c("se", "sd", 
    "conf.int", "none"), level = 0.95, xlab = deparse(substitute(factor1)), 
    ylab = paste("mean of", deparse(substitute(response))), 
	legend.lab = deparse(substitute(factor2)),  main = "",
     pch = 1:n.levs.2, lty = 1:n.levs.2, col= 1:n.levs.2, 
	levs="", levs2="", x.posn=0, y.posn=0, bty="l", l.bty="n", hatwidth = 0.3, ...)
{
    if (!is.numeric(response)) 
        stop("Argument response must be numeric.")
    xlab
    ylab
    legend.lab
    error.bars <- match.arg(error.bars)
    if (missing(factor2)) {
        if (!is.factor(factor1)) 
            stop("Argument factor1 must be a factor.")
        valid <- complete.cases(factor1, response)
        factor1 <- factor1[valid]
        response <- response[valid]
        means <- tapply(response, factor1, mean)
        sds <- tapply(response, factor1, sd)
        ns <- tapply(response, factor1, length)
        if (error.bars == "se") 
            sds <- sds/sqrt(ns)
        if (error.bars == "conf.int") 
            sds <- qt((1 - level)/2, df = ns - 1, lower.tail = FALSE) * 
                sds/sqrt(ns)
        yrange <- if (error.bars != "none") 
            c(0.9*min(means - sds),1.1*max(means + sds))
        else range(means)
	  if (levs[1] == "")
	 levs <- levels(factor1)
        n.levs <- length(levels(factor1))
        plot(c(1, n.levs), xlim=c(0.5,n.levs+.5),yrange, type = "n", 
            xlab = xlab, ylab = ylab, 
            axes = FALSE, main = main)
        points(1:n.levs, means, type = "p")
        box(bty=bty)
        axis(2)
        axis(1, at = 1:n.levs, labels = levs)
        if (error.bars != "none") 
            arrows(1:n.levs, means - sds, 1:n.levs, means + sds, 
                angle = 90, lty = 1, code = 3, length = hatwidth)
    }
    else {
        if (!(is.factor(factor1) | is.factor(factor2))) 
            stop("Arguments factor1 and factor2 must be factors.")
        valid <- complete.cases(factor1, factor2, response)
        factor1 <- factor1[valid]
        factor2 <- factor2[valid]
        response <- response[valid]
        means <- tapply(response, list(factor1, factor2), mean)
        sds <- tapply(response, list(factor1, factor2), sd)
        ns <- tapply(response, list(factor1, factor2), length)
        if (error.bars == "se") 
            sds <- sds/sqrt(ns)
        if (error.bars == "conf.int") 
            sds <- qt((1 - level)/2, df = ns - 1, lower.tail = FALSE) * 
                sds/sqrt(ns)
        yrange <- if (error.bars != "none") 
            c(0.9*min(means - sds), 1.1*max(means + sds))
        else range(means)
	  if (levs[1] == "")
		levs=levels(factor1)
        levs.1 <- levels(factor1)
        levs.2 <- levels(factor2)
        n.levs.1 <- length(levs.1)
        n.levs.2 <- length(levs.2)
          plot(c(1, n.levs.1 + .5), yrange, type = "n", xlab = xlab, 
            ylab = ylab, axes = FALSE, main = main)
        box(bty=bty)
        axis(2)
        axis(1, at = 1:n.levs.1+n.levs.2*.05, labels = levs)
       
        for (i in 1:n.levs.2) {
            points(1:n.levs.1+0.1*(i-1), means[, i], type = "p", pch = pch[i], 
             col = col[i],lty = lty[i])
            if (error.bars != "none") 
                arrows(1:n.levs.1+.1*(i-1), means[, i] - sds[, i], 
			1:n.levs.1+.1*(i-1), col = col[i],
                  means[, i] + sds[, i], angle = 90, code = 3, 
                   lty = lty[i], length = hatwidth)
        }
	  if (x.posn==0)
        x.posn <- n.levs.1 + 0.3
	  if (y.posn==0)
        y.posn <- sum(c(0.1, 0.9) * par("usr")[c(3, 4)])
        text(x.posn, y.posn, legend.lab, adj = c(0, -0.5))
	  if (levs2[1]=="") levs2=levs.2
        legend(x.posn, y.posn, levs2, pch = pch, col=col,
            lty = lty, bty=l.bty)
   }
    invisible(NULL)
} 

 

#Matched plot-------- 
# where the response is plotted with connected lines by grouping variable

matched <- function(formula, data, subset = NULL, space1 = 0.25, space2 = 0.25,  label = TRUE, xlim = NULL, ylab = NULL, xlab = NULL, enlarge.overlapped = TRUE, lwd.cor = 1, ...) 
{
  if (missing(formula)) stop("'formula' missing or incorrect")
  modf <- match.call()
  if (missing(data)) data <- environment(formula)
  ## evaluate and install the model frame
  m <- match(c("formula", "data", "subset", "weights", "na.action", "offset"),
             names(modf), 0)
  modf <- modf[c(1, m)]
  modf$drop.unused.levels <- TRUE
  modf[[1]] <- as.name("model.frame")
  mf <- eval(modf, parent.frame())
  mf <- mf[order(mf[[3]], mf[[2]]),]
  resp <- model.response(mf)
  effect <- factor(mf[[2]])
  group  <-  as.numeric(effect)
  nlev <- length(levels(effect))
  levs <- levels(factor(group))
  id <- factor(mf[[3]])   
  if (label == TRUE) space2 = space2 + 0.25
{if (is.null(xlim)) xlm <- c(1-space1, nlev+space2)
 else xlm <- xlim}
{if (is.null(xlab)) xlb <- names(mf)[2]
 else xlb <- xlab}
{if (is.null(ylab)) ylb <- names(mf)[1]
 else ylb <- ylab}
  plot(resp~group, type = "n", xaxt = "n", 
       xlim = xlm, ylab = ylb, xlab = xlb, ...)   
  axis(1, at = levs, labels = levels(effect))
  D <- data.frame(id, x = group, y = resp)
  if (any(duplicated(D)))   
  {warning("ID does not define unique cases, changing the line width for overlapping segments will not work. 
           If data contains replicated cases, it may be better to use an ID that defines unique cases, i.e. 
           'ID-replicate1' 'ID-replicate2', possibly using the function 'interaction'")
   enlarge.overlapped <- FALSE
  }
  Dspl <- split(D, D$id)
  sapply(Dspl, function(x) points(x$x, x$y, ...)) 
  segments <- lapply(1: (nlev-1), function(x) data.frame(x0 = NA, x1 = NA, y0 = NA, y1 = NA))
  
  for (i  in 1:(nlev-1)) 
  {
    wide <- data.frame(x0 = sapply(Dspl, function(x) x$x[i]), x1 = sapply(Dspl, function(x) x$x[i+1]),
                       y0 = sapply(Dspl, function(x) x$y[i]), y1 = sapply(Dspl, function(x) x$y[i+1]))
    recs <- do.call("paste", c(wide, sep = "\r"))
    same <- table(recs) #indicates the number of identical records (i.e. segments to draw)
    wide$same <- as.numeric(same)[match(recs,names(same))]
    wide$same[is.na(wide$x0) | is.na(wide$x1)] <- NA
    segments[[i]] <- wide  
  }
  {if (enlarge.overlapped)
    sapply(segments, function(x) 
      segments(x$x0, x$y0, x$x1, x$y1, lwd = x$same*lwd.cor, ...))
   else
     sapply(segments, function(x) 
       segments(x$x0, x$y0, x$x1, x$y1, ...))
  }
  label.y  <- sapply(Dspl, function(x) tail(x$y, 1L))
  label.y <- ifelse(duplicated(label.y), jitter(label.y), label.y)
  if (label) text(x=(nlev + space2/10), y = label.y, labels = names(Dspl), pos = 4)
}







# Bar and line plots with se error bars-------

bars <- function(formula, data, subset = NULL, hatwidth = 0.5, ylim = NULL,
                 xlab = NULL, ylab = NULL, ...) 
{
  dots <- list(...)
  mf <- match.call()
  if (missing(data)) data <- environment(formula)
  ## evaluate and install the model frame
  m <- match(c("formula", "data", "subset", "weights", "na.action", "offset"),
             names(mf), 0)
  mf <- mf[c(1, m)]
  mf$drop.unused.levels <- TRUE
  mf[[1]] <- as.name("model.frame")
  fr <- eval(mf, parent.frame()) 
  resp <- model.response(fr)
  effect <- factor(fr[[2]])
  {if (is.null(xlab)) xlb = names(fr)[2]
   else xlb = xlab}
  {if (is.null(ylab)) ylb = names(fr)[1]
   else ylb = ylab}
  
  
  {if (dim(fr)[2] == 3) 
    { 
      eff2 <- factor(fr[[3]])
      y = tapply(resp, list(eff2, effect), mean, na.rm = T)
      sem = tapply(resp, list(eff2, effect), se)
      ymax = max(y+sem)*1.1 
      {if (is.null(ylim)) yl = c(0, ymax)
       else yl = ylim  
      } 
      bp = barplot(y, beside = T, ylim = yl, xlab = xlb, ylab = ylb, ...)
      hw = hatwidth*0.25
      arrows(bp, y+sem, bp, y-sem, angle = 90, code = 3, length = hw)
    }
   else if (dim(fr)[2] == 2) 
     {
      y = tapply(resp, effect, mean, na.rm = T)
      sem = tapply(resp, effect, se)
      ymax = max(y+sem)*1.1
      {if (is.null(ylim)) yl = c(0, ymax)
       else yl = ylim  
      } 
      bp = barplot(y, ylim = yl, xlab = xlb, ylab = ylb, ...)
      arrows(bp, y+sem, bp, y-sem, angle = 90, code = 3, 
             length = hatwidth*1)
      }
   else stop("please use formula with one response and one or two effects, e.g. response~effect + effect2")
  }
  box(bty = "l") 
}

# Lineplot now accepts factors as the first effect argument. However, when the
# first factor is factor, but has numeric values, the function uses these values instead
# of the factor levels. So for example, if time is like 1,2,4, then it will plot as 1,2,3,4
# on the x axis with no values for 3. Otherwise, it would skip 3. If this is the desired effect
# than it makes sense to rename the x variable such as something1,2,4.
lineplot <- function(formula, data, subset = NULL, hatwidth = 0.5, type = "b", pch = 19, ylim = NULL, xlab = NULL, ylab = NULL, spread = TRUE, off = 0.1, ...) 
{
  dots <- list(...)
  mf <- match.call()
  if (missing(data)) data <- environment(formula)
  ## evaluate and install the model frame
  m <- match(c("formula", "data", "subset", "weights", "na.action", "offset"),
             names(mf), 0)
  mf <- mf[c(1, m)]
  mf$drop.unused.levels <- TRUE
  mf[[1]] <- as.name("model.frame")
  fr <- eval(mf, parent.frame()) 
  resp <- model.response(fr)
  effect <- fr[[2]]
  nlev1 = length(unique(effect))
  {if (is.null(xlab)) xlb = names(fr)[2]
   else xlb = xlab}
  {if (is.null(ylab)) ylb = names(fr)[1]
   else ylb = ylab}
  {if (dim(fr)[2] == 3) 
    { 
      eff2 <- factor(fr[[3]])
      nlev2 <- length(levels(eff2))
      matr = tapply(resp, list(effect, eff2), mean, na.rm = T)
      sem = tapply(resp, list(effect, eff2), se)
    }
   else if (dim(fr)[2] == 2) {
     matr = tapply(resp, effect, mean, na.rm = T)
     sem = tapply(resp, effect, se)
     nlev2 = 1
    }                
   else stop("please use formula with one response and one or two effects, e.g. response~effect + effect2")
  }
  if (!is.numeric(type.convert(dimnames(matr)[[1]], as.is = T))){
    if (is.factor(effect)) {
      x = 1:nlevels(effect)
      names(x) = levels(effect)
      xaxt = "n"
    }
    else {
      effect = factor(effect)
      x = 1:nlevels(effect)
      names(x) = levels(effect)
      xaxt = "n"
    }
  }
  else {
    x = as.numeric(dimnames(matr)[[1]])
    xaxt = "s"
  }
  y = matr
  ymin = min(y-sem)*0.9
  ymax = max(y+sem)*1.1
  hw = hatwidth*0.25
  #check if col was specified in the main call
  if (is.null(dots$col)) color = rep(1:nlev2, each = nlev1)
  else color = rep(dots$col, each = nlev1)  
   
  #check if ylim was specified in the main call
  if (is.null(ylim)) yl = c(ymin, ymax)
  else yl = ylim  
   
  {if (spread)
    {
      dif = mean(diff(x))*off
      posit = 1:nlev2 - median(1:nlev2)
      offs = posit * dif
      X =  sapply(offs, function(o) o+x)
      matplot(x=X, y=y, type = type, ylim = yl, pch = pch, 
              xlab = xlb, ylab = ylb, xaxt = xaxt, ...)
      arrows(X, y+sem, X, y-sem, angle = 90, code = 3, length = hw, 
             col = color)
    }   
   else  
    { 
      matplot(x=x, y=y, type = type, ylim = yl, pch = pch, 
              xlab = xlb, ylab = ylb, xaxt = xaxt, ...)
      arrows(x, y+sem, x, y-sem, angle = 90, code = 3, length = hw, 
             col = color)
    }
  }
  if (xaxt == "n") {
    axis(1, at = x, labels = names(x))
  }
    
  box(bty = "l") 
}







pairs2 <- function (x,y,smooth=TRUE, digits=2, prefix="", cex.cor=NULL, stars = TRUE, ...)
{
  panel.cor <- function(x, y, ...)
  {
    usr <- par("usr"); on.exit(par(usr))
    par(usr = c(0, 1, 0, 1))
    r.obj = cor.test(x, y,use="pairwise",...)
    r = as.numeric(r.obj$estimate)
    p = r.obj$p.value
    {if (stars) mystars <- 
       ifelse(p < .001, "***", ifelse(p < .01, "** ", ifelse(p < .05, "* ", ifelse(p < 0.1, ".", " "))))
    else mystars  <-  NULL  }
    txt <- format(c(r, 0.123456789), digits=digits)[1]
    txt <- paste(prefix, txt, mystars, sep="")
    cex = ifelse(is.null(cex.cor), 0.8/strwidth(txt), cex.cor)
    cexfinal = ifelse(stars, 1.2, cex * abs(r)+0.5)
    text(0.5, 0.5, txt, cex = cexfinal)
  }
  panel.hist <- function(x)
  {
    usr <- par("usr"); on.exit(par(usr))
    par(usr = c(usr[1:2], 0, 1.5) )
    h <- hist(x, plot = FALSE)
    breaks <- h$breaks; nB <- length(breaks)
    y <- h$counts; y <- y/max(y)
    rect(breaks[-nB], 0, breaks[-1], y, col="cyan")
  }
  {if (smooth)
    pairs(x,diag.panel=panel.hist,lower.panel=panel.cor,upper.panel=panel.smooth)
  else #smooth is not true
    pairs(x,diag.panel=panel.hist,lower.panel=panel.cor)}
} 


# Table formatting -----
fixed_digits <- function(xs, n = 2) {
  formatC(xs, digits = n, format = "f")
}

# nice p-values for tables...
nice_pval <- function(ps, html = FALSE) {
  tiny <- ifelse(html, "&lt;&nbsp;0.001", "<0.001")
  ps_chr <- fixed_digits(ps, 3) 
  ps_chr[ps < 0.001] <- tiny
  ps_chr
}

# ...and for in-text
nice_pval2 <- function(ps, html = FALSE) {
  tiny <- ifelse(html, "&lt;&nbsp;0.001", "<0.001")
  ps_chr <- paste0("= ", fixed_digits(ps, 3))
  ps_chr[ps < 0.001] <- tiny
  ps_chr
}


nice_modtable <- function(table, stat.digit = 2, pval = "p.value"){
  if (is.matrix(table)) table <-  as.data.frame(table)
  names(table)[1] <-  "term"
  col.d2 <-  which(!names(table) %in% c("term", pval))
  for (i in col.d2){
    table[,i] <-  fixed_digits(table[,i], stat.digit)
  }
  table[,names(table) %in% pval] <- nice_pval(table[,names(table) %in% pval])
  table
}

nicetab <- function(x) knitr::kable(nice_modtable(broom::tidy(x)))

