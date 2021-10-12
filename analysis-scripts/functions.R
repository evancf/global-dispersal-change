# Various functions used in the analysis and making custom plots.

ipak <- function(pkg){
  new.pkg <- pkg[!(pkg %in% installed.packages()[, "Package"])]
  if (length(new.pkg))
    install.packages(new.pkg, dependencies = TRUE)
  sapply(pkg, require, character.only = TRUE)
}

l.unique <- function(x) length(unique(x))

as.num.char <- function(x) as.numeric(as.character(x))

# A couple functions that help go from long format to bipartite matrix format
net.spread.inside <- function(split.by = split.by, split.vals = split.vals, 
                              tax.type = tax.type, data.type = data.type,
                              long.df = long.df){
  
  x <- long.df[which(long.df[,split.by] == split.vals),]
  a <- x[ ,paste("animal", tax.type, sep = ".")]
  p <- x[ ,paste("plant", tax.type, sep = ".")]
  identifier <- paste(a, p, sep = " ~ ")
  y <- tapply(x$value, identifier, sum)
  
  if(data.type == "bin"){
    y <- ifelse(y > 0, 1, 0)
  }
  
  z <- spread(as.data.frame(cbind(str_split_fixed(names(y), " ~ ", 2), y)),
              value = 3,
              key = 1,
              fill = 0)
  #print(split.vals) # In case you want to track down problems
  
  rownames(z) <- z[,1]
  
  z <- z[which(!rownames(z) == "NA"), ]
  z <- z[, which(!colnames(z) == "NA")]
  
  
  # This chunk was added to make it so cases where there is 1 or 0 of the 
  # focal plant/frug taxa (at whatever focal taxa level, e.g., accepted species)
  # the first column doesn't have to get removed (because it's sort of irrelevant because these will not be added to the list).
  if(!is.null(dim(z))){
    
    z <- z[,-1] # Can do z[,-1, drop = F] to retain a 1 column dataframe rather than vector
    
    if(!is.null(dim(z))){ # The next chunk wont work if there is only 1 
      z[] <- mutate_all(z, list(as.num.char)) # Make numeric while retaining row names...
      
      z <- z[which(!rowSums(z) == 0), ] # Get rid of plant or animal taxa with no interactions
      z <- z[, which(!colSums(z) == 0)]
    }
  }
  
  z
  
}

net.spread <- function(split.by, split.vals, 
                       tax.type, data.type, 
                       long.df = long.df,
                       min.taxa = 2){
  if(min.taxa <= 1) stop("cannot make network with only one plant or animal taxon") # This may be obvious
  zz <- apply(cbind(split.by, as.character(split.vals), tax.type, data.type), 1, 
              function(xx) net.spread.inside(split.by = xx[1], split.vals = xx[2], tax.type = xx[3], data.type = xx[4], long.df = long.df))
  names(zz) <- split.vals
  zz <- zz[!unlist(lapply(zz, function(zzz) any(dim(zzz) <= min.taxa) | is.null(dim(zzz))))]
  zz
}





# Here's a slight change to fourfoldplot

my.fourfoldplot <- function (x, color = c("#99CCFF", "#6699CC"), conf.level = 0.95, 
                             std = c("margins", "ind.max", "all.max"), margin = c(1, 2), 
                             space = 0.2, main = NULL, mfrow = NULL, mfcol = NULL,
                             margin.text = T, cex.num = 1) 
{
  if (!is.array(x)) 
    stop("'x' must be an array")
  if (length(dim(x)) == 2L) {
    x <- if (is.null(dimnames(x))) 
      array(x, c(dim(x), 1L))
    else array(x, c(dim(x), 1L), c(dimnames(x), list(NULL)))
  }
  if (length(dim(x)) != 3L) 
    stop("'x' must be 2- or 3-dimensional")
  if (any(dim(x)[1L:2L] != 2L)) 
    stop("table for each stratum must be 2 by 2")
  dnx <- dimnames(x)
  if (is.null(dnx)) 
    dnx <- vector("list", 3L)
  for (i in which(sapply(dnx, is.null))) dnx[[i]] <- LETTERS[seq_len(dim(x)[i])]
  if (is.null(names(dnx))) 
    i <- 1L:3L
  else i <- which(is.null(names(dnx)))
  if (any(i)) 
    names(dnx)[i] <- c("Row", "Col", "Strata")[i]
  dimnames(x) <- dnx
  k <- dim(x)[3L]
  if (!((length(conf.level) == 1) && is.finite(conf.level) && 
        (conf.level >= 0) && (conf.level < 1))) 
    stop("'conf.level' must be a single number between 0 and 1")
  if (conf.level == 0) 
    conf.level <- FALSE
  std <- match.arg(std)
  findTableWithOAM <- function(or, tab) {
    m <- rowSums(tab)[1L]
    n <- rowSums(tab)[2L]
    t <- colSums(tab)[1L]
    if (or == 1) 
      x <- t * n/(m + n)
    else if (or == Inf) 
      x <- max(0, t - m)
    else {
      A <- or - 1
      B <- or * (m - t) + (n + t)
      C <- -t * n
      x <- (-B + sqrt(B^2 - 4 * A * C))/(2 * A)
    }
    matrix(c(t - x, x, m - t + x, n - x), nrow = 2)
  }
  drawPie <- function(r, from, to, n = 500, col = NA) {
    p <- 2 * pi * seq.int(from, to, length.out = n)/360
    x <- c(cos(p), 0) * r
    y <- c(sin(p), 0) * r
    polygon(x, y, col = col)
    invisible(NULL)
  }
  stdize <- function(tab, std, x) {
    if (std == "margins") {
      if (all(sort(margin) == c(1L, 2L))) {
        u <- sqrt(odds(tab)$or)
        u <- u/(1 + u)
        y <- matrix(c(u, 1 - u, 1 - u, u), nrow = 2L)
      }
      else if (margin %in% c(1, 2)) 
        y <- prop.table(tab, margin)
      else stop("incorrect 'margin' specification")
    }
    else if (std == "ind.max") 
      y <- tab/max(tab)
    else if (std == "all.max") 
      y <- tab/max(x)
    y
  }
  odds <- function(x) {
    if (length(dim(x)) == 2L) {
      dim(x) <- c(dim(x), 1L)
      k <- 1
    }
    else k <- dim(x)[3L]
    or <- double(k)
    se <- double(k)
    for (i in 1:k) {
      f <- x[, , i]
      storage.mode(f) <- "double"
      if (any(f == 0)) 
        f <- f + 0.5
      or[i] <- (f[1L, 1L] * f[2L, 2L])/(f[1L, 2L] * f[2L, 
                                                      1L])
      se[i] <- sqrt(sum(1/f))
    }
    list(or = or, se = se)
  }
  gamma <- 1.25
  debug <- FALSE
  angle.f <- c(90, 180, 0, 270)
  angle.t <- c(180, 270, 90, 360)
  opar <- par(mar = c(0, 0, if (is.null(main)) 0 else 2.5, 
                      0))
  on.exit(par(opar))
  byrow <- FALSE
  if (!is.null(mfrow)) {
    nr <- mfrow[1L]
    nc <- mfrow[2L]
  }
  else if (!is.null(mfcol)) {
    nr <- mfcol[1L]
    nc <- mfcol[2L]
    byrow <- TRUE
  }
  else {
    nr <- ceiling(sqrt(k))
    nc <- ceiling(k/nr)
  }
  if (nr * nc < k) 
    stop("incorrect geometry specification")
  if (byrow) 
    indexMatrix <- expand.grid(1:nc, 1:nr)[, c(2, 1)]
  else indexMatrix <- expand.grid(1:nr, 1:nc)
  totalWidth <- nc * 2 * (1 + space) + (nc - 1L) * space
  totalHeight <- if (k == 1) 
    2 * (1 + space)
  else nr * (2 + (2 + gamma) * space) + (nr - 1L) * space
  xlim <- c(0, totalWidth)
  ylim <- c(0, totalHeight)
  dev.hold()
  on.exit(dev.flush(), add = TRUE)
  plot.new()
  plot.window(xlim = xlim, ylim = ylim, asp = 1)
  o <- odds(x)
  scale <- space/(2 * strheight("Ag"))
  v <- 0.95 - max(strwidth(as.character(c(x)), cex = scale))/2
  for (i in 1:k) {
    tab <- x[, , i]
    fit <- stdize(tab, std, x)
    xInd <- indexMatrix[i, 2L]
    xOrig <- 2 * xInd - 1 + (3 * xInd - 2) * space
    yInd <- indexMatrix[i, 1L]
    yOrig <- if (k == 1) 
      (1 + space)
    else (totalHeight - (2 * yInd - 1 + ((3 + gamma) * yInd - 
                                           2) * space))
    plot.window(xlim - xOrig, ylim - yOrig, asp = 1)
    if (debug) {
      abline(h = -1 - space)
      abline(h = 1 + space)
      abline(h = 1 + (1 + gamma) * space)
      abline(v = -1 - space)
      abline(v = 1 + space)
    }
    u <- 1 + space/2
    adjCorr <- 0.2
    if(margin.text){
      text(0, u, paste(names(dimnames(x))[1L], dimnames(x)[[1L]][1L], 
                       sep = ": "), adj = c(0.5, 0.5 - adjCorr), cex = scale)
      text(-u, 0, paste(names(dimnames(x))[2L], dimnames(x)[[2L]][1L], 
                        sep = ": "), adj = c(0.5, 0.5 - adjCorr), cex = scale, 
           srt = 90)
      text(0, -u, paste(names(dimnames(x))[1L], dimnames(x)[[1L]][2L], 
                        sep = ": "), adj = c(0.5, 0.5 + adjCorr), cex = scale)
      text(u, 0, paste(names(dimnames(x))[2L], dimnames(x)[[2L]][2L], 
                       sep = ": "), adj = c(0.5, 0.5 + adjCorr), cex = scale, 
           srt = 90)
    }
    
    if (k > 1) {
      text(0, 1 + (1 + gamma/2) * space, paste(names(dimnames(x))[3L], 
                                               dimnames(x)[[3L]][i], sep = ": "), cex = gamma * 
             scale)
    }
    d <- odds(tab)$or
    drawPie(sqrt(fit[1, 1]), 90, 180, col = color[1 + (d > 
                                                         1)])
    drawPie(sqrt(fit[2, 1]), 180, 270, col = color[2 - (d > 
                                                          1)])
    drawPie(sqrt(fit[1, 2]), 0, 90, col = color[2 - (d > 
                                                       1)])
    drawPie(sqrt(fit[2, 2]), 270, 360, col = color[1 + (d > 
                                                          1)])
    u <- 1 - space/2
    text(c(-v, -v, v, v), c(u, -u, u, -u), as.character(c(tab)), 
         cex = scale*cex.num)
    if (is.numeric(conf.level)) {
      or <- o$or[i]
      se <- o$se[i]
      theta <- or * exp(stats::qnorm((1 - conf.level)/2) * 
                          se)
      tau <- findTableWithOAM(theta, tab)
      r <- sqrt(c(stdize(tau, std, x)))
      for (j in 1:4) drawPie(r[j], angle.f[j], angle.t[j])
      theta <- or * exp(stats::qnorm((1 + conf.level)/2) * 
                          se)
      tau <- findTableWithOAM(theta, tab)
      r <- sqrt(c(stdize(tau, std, x)))
      for (j in 1:4) drawPie(r[j], angle.f[j], angle.t[j])
    }
    #polygon(c(-1, 1, 1, -1), c(-1, -1, 1, 1))
    lines(c(-1, 1), c(0, 0))
    for (j in seq.int(from = -0.8, to = 0.8, by = 0.2)) lines(c(j, 
                                                                j), c(-0.02, 0.02))
    for (j in seq.int(from = -0.9, to = 0.9, by = 0.2)) lines(c(j, 
                                                                j), c(-0.01, 0.01))
    lines(c(0, 0), c(-1, 1))
    for (j in seq.int(from = -0.8, to = 0.8, by = 0.2)) lines(c(-0.02, 
                                                                0.02), c(j, j))
    for (j in seq.int(from = -0.9, to = 0.9, by = 0.2)) lines(c(-0.01, 
                                                                0.01), c(j, j))
  }
  if (!is.null(main)) 
    mtext(main, cex = 1.5, adj = 0.5)
  return(invisible())
}



# Coppied from here: https://logfc.wordpress.com/2017/03/15/adding-figure-labels-a-b-c-in-the-top-left-corner-of-the-plotting-region/

fig_label <- function(text, cex = 1, font = 1){
  pu <- par("usr")
  text(text, pu[1], pu[4], cex = cex, font = font, xpd = T)
}

put.fig.letter <- function(label, location="topleft", x=NULL, y=NULL, 
                           offset=c(0, 0), ...) {
  if(length(label) > 1) {
    warning("length(label) > 1, using label[1]")
  }
  if(is.null(x) | is.null(y)) {
    coords <- switch(location,
                     topleft = c(0.15,0.9),
                     topcenter = c(0.5525,0.9),
                     # topright = c(0.985, 0.98),
                     # bottomleft = c(0.015, 0.02), 
                     # bottomcenter = c(0.5525, 0.02), 
                     # bottomright = c(0.985, 0.02),
                     c(0.15, 0.9) )
  } else {
    coords <- c(x,y)
  }
  this.x <- grconvertX(coords[1] + offset[1], from="nfc", to="user")
  this.y <- grconvertY(coords[2] + offset[2], from="nfc", to="user")
  text(labels=label[1], x=this.x, y=this.y, xpd=T, ...)
}




addImg <- function(
  obj, # an image file imported as an array (e.g. png::readPNG, jpeg::readJPEG)
  x = NULL, # mid x coordinate for image
  y = NULL, # mid y coordinate for image
  width = NULL, # width of image (in x coordinate units)
  interpolate = TRUE # (passed to graphics::rasterImage) A logical vector (or scalar) indicating whether to apply linear interpolation to the image when drawing. 
){
  if(is.null(x) | is.null(y) | is.null(width)){stop("Must provide args 'x', 'y', and 'width'")}
  USR <- par()$usr # A vector of the form c(x1, x2, y1, y2) giving the extremes of the user coordinates of the plotting region
  PIN <- par()$pin # The current plot dimensions, (width, height), in inches
  DIM <- dim(obj) # number of x-y pixels for the image
  ARp <- DIM[1]/DIM[2] # pixel aspect ratio (y/x)
  WIDi <- width/(USR[2]-USR[1])*PIN[1] # convert width units to inches
  HEIi <- WIDi * ARp # height in inches
  HEIu <- HEIi/PIN[2]*(USR[4]-USR[3]) # height in units
  rasterImage(image = obj, 
              xleft = x-(width/2), xright = x+(width/2),
              ybottom = y-(HEIu/2), ytop = y+(HEIu/2), 
              interpolate = interpolate)
}




pointerfun <- function(spp, pic.x, pic.y, response.var = "dist", log.y = T,
                       col = "grey",
                       single.plant = NA){
  
  if(is.na(single.plant)){
    spp.xy <- apply(subset(fig1.set, 
                           animal.phylo.id == spp)[,c("body_mass_median",
                                                      response.var)], 2, mean)
  } else{
    spp.xy <- apply(subset(fig1.set, 
                           animal.phylo.id == spp &
                             plant.phylo.id == single.plant)[,c("body_mass_median",
                                                                response.var)], 2, mean)
  }
  
  
  if(log.y == T){
    spp.xy[2] <- log(spp.xy[2])
  }
  segments(spp.xy[1], spp.xy[2],
           pic.x, pic.y,
           col = col,
           lwd = pic.lwd)
  draw.circle(spp.xy[1], spp.xy[2], 
              radius = pointer.radius,
              lwd = pic.lwd * 0.66,
              col = "grey40",
              border = col)
  draw.circle(pic.x, pic.y,
              radius = pic.radius,
              lwd = pic.lwd*2,
              border = col)
}






net.plot <- function(net, grey.breaks = NA, 
                     show.net = NA, border = NA,
                     block.out.introd = F,
                     introd.sp.net.to.add = NULL,
                     introd.spp = NA,
                     introd.breaks = NA,
                     sp.cex =  0.7, col.text = T, row.text = T,
                     obs.rect.net = NA, 
                     inner.rect.col = "lightblue",
                     xlim.max = NA,
                     ylim.max = NA,
                     supress.legend = F,
                     leg.title = NULL, exp.leg = F,
                     leg.x.pos = 0.7,
                     leg.y.pos = 0.63,  # 0.69
                     leg.cex = 0.8,
                     sig.digits = 2, 
                     sp.abbr = F,
                     col.srt = 45, ...){
  
  if(!is.na(show.net)){
    net <- net * ((round(show.net * grey.breaks) / grey.breaks) != 0)
  } 
  
  inds <- which(net != 0, arr.ind = T)
  
  col.vals <- rgb(0,0,0,
                  unlist(net[inds] > 0))
  
  if(!is.na(grey.breaks)){
    rounded.vals <- round(unlist(net[inds])/max(unlist(net[inds])) * grey.breaks) / grey.breaks
    col.vals <- rgb(0,0,0,
                    rounded.vals)
  }
  
  
  plot(NA,
       asp = 1,
       xlim = c(1, ifelse(!is.na(xlim.max), xlim.max, max(inds[,2]) + 1)),
       ylim = -c(ifelse(!is.na(ylim.max), ylim.max, max(inds[,1]) + 1), 1),
       frame = F,
       xaxt = "n", yaxt = "n",
       xlab = "", ylab = "",...)
  
  rect(inds[,2], -inds[,1],  inds[,2]+1, -inds[,1]-1,
       col = col.vals,
       border = border)
  
  
  if(block.out.introd){
    
    if(any(colnames(net) %in% introd.spp)){
      rect(which(colnames(net) %in% introd.spp), -1, 
           which(colnames(net) %in% introd.spp) +  1, -max(inds[,1]) - 1,
           col = "lightblue", border = border)
    }
    
    if(any(rownames(net) %in% introd.spp)){
      rect(1, -which(rownames(net) %in% introd.spp),  
           max(inds[,2]) + 1,  -which(rownames(net) %in% introd.spp) -  1,
           col = "lightblue", border = border)
    }
    
  }
  
  if(!is.null(introd.sp.net.to.add)){
    
    if(any(colnames(net) %in% introd.spp)){
      rect(which(colnames(net) %in% introd.spp), -1, 
           which(colnames(net) %in% introd.spp) +  1, -max(inds[,1]) - 1,
           col = "white", border = border)
    }
    
    if(any(rownames(net) %in% introd.spp)){
      rect(1, -which(rownames(net) %in% introd.spp),  
           max(inds[,2]) + 1,  -which(rownames(net) %in% introd.spp) -  1,
           col = "white", border = border)
    }
    introd.net <- introd.sp.net.to.add * ifelse(apply(expand.grid(rownames(introd.sp.net.to.add) %in% introd.spp, 
                                                                  colnames(introd.sp.net.to.add) %in% introd.spp), 1, 
                                                      function(x) x[1] | x[2]) %>% matrix(nrow=nrow(net)),
                                                1, 0)
    
    introd.inds <- which(introd.net != 0, arr.ind = T)
    
    if(!is.na(introd.breaks)){
      introd.rounded.vals <- round(unlist(introd.net[introd.inds])/max(unlist(introd.net[introd.inds])) * introd.breaks) / introd.breaks
      introd.col.vals <- rgb(0,0,0.9,
                             introd.rounded.vals)
    }
    
    rect(introd.inds[,2], -introd.inds[,1],  introd.inds[,2]+1, -introd.inds[,1]-1,
         col = introd.col.vals,
         border = border)
    
  }
  
  if(!is.na(obs.rect.net)){
    obs.inds <- which(obs.rect.net != 0, arr.ind = T)
    rect(obs.inds[,2]+0.3, -obs.inds[,1]-0.3,  obs.inds[,2]+0.7, -obs.inds[,1]-0.7,
         col = inner.rect.col,
         border = NA,
         lwd = 0.5)
  }
  
  if(col.text){
    if(sp.abbr == T){
      col.sp.names <- paste(substr(word(colnames(net), 1), 1, 4), 
                            ". ", 
                            substr(word(colnames(net), 2), 1, 2),
                            ".",
                            sep = "")
    } else{
      col.sp.names <- paste(word(colnames(net), 1), 
                            " ", 
                            substr(word(colnames(net), 2), 1, 1),
                            ".",
                            sep = "")
    }
    text(1:max(inds[,2])-0.5, -0.5, col.sp.names, cex = sp.cex,
         col = ifelse(colnames(net) %in% introd.spp, "blue", "black"),
         pos = 4,
         font = 3,
         xpd = T,
         srt = col.srt)
  }
  
  if(row.text){
    if(sp.abbr == T){
      row.sp.names <- paste(substr(word(rownames(net), 1), 1, 4), 
                            ". ", 
                            substr(word(rownames(net), 2), 1, 2),
                            ".",
                            sep = "")
    } else{
      row.sp.names <- paste(word(rownames(net), 1), 
                            " ", 
                            substr(word(rownames(net), 2), 1, 1),
                            ".",
                            sep = "")
    }
    text(1, -0.5-(1:max(inds[,1])), row.sp.names, cex = sp.cex,
         col = ifelse(rownames(net) %in% introd.spp, "blue", "black"),
         pos = 2,
         font = 3,
         xpd = T)
  }
  
  
  if(supress.legend == F){
    
    if(!is.na(grey.breaks)){
      leg.vals <- signif(seq(from = 0, 
                             to = max(unlist(net[inds])), 
                             length.out = grey.breaks), sig.digits)[-1]
      if(exp.leg){
        leg.labels <- signif(exp(seq(from = 0, 
                                     to = max(unlist(net[inds])), 
                                     length.out = grey.breaks)), sig.digits)[-1]
      } else{
        leg.labels <- leg.vals
      }
      
      legend(max(inds[,2])*leg.x.pos,
             -max(inds[,1])*leg.y.pos,
             legend = leg.labels,
             title = leg.title,
             cex = leg.cex,
             title.adj = 0,
             bty = "n",
             pch = 15,
             pt.cex = 1,
             col = rgb(0,0,0, leg.vals / max(leg.vals)))
    }
    
    if(!is.na(introd.breaks)){
      leg.vals <- signif(seq(from = 0, 
                             to = max(unlist(introd.net[introd.inds])), 
                             length.out = introd.breaks), sig.digits)[-1]
      if(exp.leg){
        leg.labels <- signif(exp(seq(from = 0, 
                                     to = max(unlist(introd.net[introd.inds])), 
                                     length.out = introd.breaks)), sig.digits)[-1]
      } else{
        leg.labels <- leg.vals
      }
      
      legend(max(introd.inds[,2])*leg.x.pos,
             -max(introd.inds[,1])*leg.y.pos,
             legend = leg.labels,
             title = leg.title,
             title.adj = 0,
             bty = "n",
             pch = 15,
             pt.cex = 1.3,
             col = rgb(0,0,1, leg.vals / max(leg.vals)))
    }
  }
  
}






net.web <- function(net, grey.breaks = NA, show.net = NA, border = NA,
                    sp.cex =  0.7, animal.text = T, plant.text = T,
                    obs.rect.net = NA, leg.title = NULL, exp.leg = F){
  
  if(!is.na(show.net)){
    net <- net * ((round(show.net * grey.breaks) / grey.breaks) != 0)
  } 
  
  inds <- which(net != 0, arr.ind = T)
  
  col.vals <- rgb(0,0,0,
                  unlist(net[inds] > 0))
  
  if(!is.na(grey.breaks)){
    rounded.vals <- round(unlist(net[inds])/max(unlist(net[inds])) * grey.breaks) / grey.breaks
    col.vals <- rgb(0,0,0,
                    rounded.vals)
  }
  
  x.wide <- max(inds[,1]) * 0.2
  
  plot(NA,
       asp = 1,
       xlim = c(0, x.wide),
       ylim = -c(max(inds[,1]) + 1, 1),
       frame = F,
       xaxt = "n", yaxt = "n",
       xlab = "", ylab = "")
  
  segments(1, -inds[,1]-0.5, x.wide, -inds[,2]-0.5,
           lwd = 2,
           col = col.vals)
  
  
  if(animal.text){
    text(x.wide, -0.5-(1:max(inds[,2])), colnames(net), cex = sp.cex,
         pos = 4,
         font = 3,
         xpd = T)
  }
  
  if(plant.text){
    text(1, -0.5-(1:max(inds[,1])), rownames(net), cex = sp.cex,
         pos = 2,
         font = 3,
         xpd = T)
  }
  
}


plot.focal.net <- function(focal.net.id, panel.lab = NA, ...){
  
  focal.set <- subset(int.set, net.id == focal.net.id)
  
  introd.spp <- unique(c(subset(focal.set, animal.native.status == "introduced")$animal.phylo.id,
                         subset(focal.set, plant.native.status == "introduced")$plant.phylo.id))
  
  
  
  focal.introd.pred.set <- focal.all.pred.set <- focal.set
  
  focal.all.pred.set$value <- predict(object = int.mod,
                                      newdata = focal.all.pred.set,
                                      n.trees = gbm.perf(int.mod, method = "cv", plot.it = F),
                                      type = "response")
  
  focal.introd.pred.set$value <- predict(object = introd.mod,
                                         newdata = focal.introd.pred.set,
                                         n.trees = gbm.perf(introd.mod, method = "cv", plot.it = F),
                                         type = "response")
  
  
  
  focal.net <- net.spread(split.by = "net.id", 
                          split.vals = focal.net.id,
                          tax.type = "phylo.id", 
                          data.type = "quant",
                          long.df = focal.set)[[1]]
  focal.net <- order.net(focal.net, focal.net)
  
  
  focal.all.pred.net <- net.spread(split.by = "net.id", 
                                   split.vals = focal.net.id,
                                   tax.type = "phylo.id", 
                                   data.type = "quant",
                                   long.df = focal.all.pred.set)[[1]]
  focal.all.pred.net <- order.net(focal.all.pred.net, focal.net)
  
  
  focal.introd.pred.net <- net.spread(split.by = "net.id", 
                                      split.vals = focal.net.id,
                                      tax.type = "phylo.id", 
                                      data.type = "quant",
                                      long.df = focal.introd.pred.set)[[1]]
  focal.introd.pred.net <- order.net(focal.introd.pred.net, focal.net)
  
  
  net.plot(focal.net,
           introd.breaks = 6,
           obs.rect.net = focal.net,
           inner.rect.col = "black",
           introd.sp.net.to.add = focal.introd.pred.net, 
           introd.spp = introd.spp,
           leg.x.pos = 0.9,
           leg.title = "Probability",
           supress.legend = T,
           ...)
  text(-1, 1, panel.lab, font = 2, cex = 1.2, xpd = T)
}





blank.map <- function(cex = 1.2, meta = metanet,
                      add.box = T,
                      map.col = rgb(0,0,0,0.1),
                      ylim = c(-60, 80), col = rgb(117,112,179, maxColorValue = 255, 200)){
  map('world',col = map.col, 
      fill = TRUE, #bg="white", 
      lwd = 0.01, border = NA,#border = rgb(0,0,0,0.1), 
      ylim = ylim,
      xlim = c(-180,180), 
      mar = c(0,0,0,0))
  if(add.box){
    box(which = "figure")
  }
  points(latitude ~ longitude, data = meta,
         col = col, cex = cex, pch = 16)
}


