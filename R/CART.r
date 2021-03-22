# ****************************************************************************
#
# Projekt:  ModTools
#
# Zweck:    some CART stuff, routines for handling rparts
#
# Autor:    partly by B. Compton, compiled and extended by Andri Signorell
#
# ****************************************************************************



# plot.rpart <- function (x, uniform = FALSE, branch = 1, compress = FALSE, nspace,
#           margin = 0, minbranch = 0.3, ...) {
#   rpart.plot(x=x, uniform = uniform, branch = branch, compress = compress, nspace = nspace,
#              margin =margin, minbranch = minbranch, ...)
# }

#
# plot.rpart <- function (x = stop("no 'x' arg"), type = 0, extra = 0, under = FALSE,
#           clip.right.labs = TRUE, fallen.leaves = FALSE, branch = if (fallen.leaves) 1 else 0.2,
#           uniform = TRUE, digits = 2, varlen = -8, faclen = 3, cex = NULL,
#           tweak = 1, compress = TRUE, ycompress = uniform, snip = FALSE,
#           ...){
#   rpart.plot (x , type , extra , under ,
#               clip.right.labs , fallen.leaves , branch ,
#               uniform , digits , varlen , faclen , cex,
#               tweak , compress , ycompress , snip,
#               ...)
# }
#



.myStripStyle <- function(which.panel, factor.levels, ...) {
  panel.rect(0, -0.5, 1, 1,
             col = "grey",
             border = 1)
  panel.text(x = 0.5, y = 0.25,
             font=2,
             lab = gettextf("%s", factor.levels[which.panel]),
             col = "black")
}





# Rules ------------------
Rules <- function(x, node=NULL, leafonly = FALSE) {

  if (!inherits(x, "rpart")) stop("Not a legitimate rpart tree")
  #
  # Get some information.
  #


  if(is.null(node))
    node <- rownames(x$frame)
  else
    node <- as.character(node)

  frm <- x$frame[node, ]

  if(leafonly)
    # the leafs
    frm     <- frm[frm$var == "<leaf>", ]

  names   <- row.names(frm)
  ylevels <- attr(x, "ylevels")
  ds.size <- x$frame[1,]$n

  if(nrow(frm) > 0)
    structure(list(
      frame=frm,
      ylevels=ylevels,
      ds.size=ds.size,
      path=rpart::path.rpart(x, nodes=as.numeric(row.names(frm)), print.it=FALSE)),
      class="Rules")
  else
    NA

}


print.Rules <- function(x, ...){

  # Print each leaf node as a rule.

  for (i in 1:nrow(x$frame)) {

    cat("\n")
    cat(sprintf(" Rule number: %s ", rownames(x$frame)[i]))

    if (x$frame[i, 1] == "<leaf>") {

      cat(sprintf("[yval=%s cover=%d (%.0f%%) prob=%0.2f]",
                  x$ylevels[x$frame[i,]$yval], x$frame[i,]$n,
                  round(100 * x$frame[i,]$n / x$ds.size), x$frame[i,]$yval2[,5]))
    }

    cat("\n")
    cat(sprintf("   %s\n", x$path[[i]][-1]), sep="")

  }

  cat("\n")

}




# TreeDepth ------------------

TreeDepth <- function (nodes)  {
    depth <- floor(log(nodes, base = 2) + 1e-7)
    as.vector(depth - min(depth))
}



# Desc.rpart ------------------

Desc.rpart <- function(x, nodes = NULL, cp = 0, digits = getOption("digits"), ...) {

  #  SCCS  @(#)summary.rpart.s  1.18 07/05/01
  #    adapted for single node summary:   22.12.2009 / Andri Signorell

  #  purpose:		prints the summary for single nodes, where nodes is the node id
  #               then rest remains original code
  #  example: 	node.summary( r.rpart, nodes=c(32) )

  object <- x

  if(!inherits(object, "rpart")) stop("Not legitimate rpart object")

  # If this is an older-style rpart object, convert it
  #  either way, rename it to "x" to save typing
  # if (!is.null(object$frame$splits)) x <- rpconvert(object)
  # else  x <- object
  # get rid of old hats
  x <- object

  # omit <- x$na.action
  # n <- x$frame$n
  # if (length(omit))
  # cat("  n=", n[1], " (", naprint(omit), ")\n\n", sep="")
  # else cat("  n=", n[1], "\n\n")

  #    print(x$cptable, digits=digits)

  ff <- x$frame
  ylevel <- attr(x,'ylevels')
  id <- as.integer(row.names(ff))
  parent.id <- ifelse(id==1,1, floor(id/2))
  parent.cp <- ff$complexity[match(parent.id, id)]
  rows <- (1:length(id))[parent.cp > cp]
  if (length(rows)>0) rows <- rows[order(id[rows])]
  else rows <- 1
  is.leaf <- (ff$var=='<leaf>')
  index <- cumsum(c(1, ff$ncompete + ff$nsurrogate + 1*(!is.leaf)))

  if(!all(is.leaf)) {  #skip these lines for a "no splits" tree
    sname <- dimnames(x$splits)[[1]]
    cuts <- vector(mode='character', length=nrow(x$splits))
    temp <- x$splits[ ,2]
    for (i in 1:length(cuts)) {
      if (temp[i] == -1)
        cuts[i] <-paste("<", format(signif(x$splits[i,4], digits=digits)))
      else if (temp[i] ==1)
        cuts[i] <-paste("<", format(signif(x$splits[i,4], digits=digits)))
      else cuts[i]<- paste("splits as ",
                           paste(c("L", "-", "R")[x$csplit[x$splits[i,4], 1:temp[i]]],
                                 collapse='', sep=''), collapse='')
    }
    # S-PLUS 4.0 can't handle null vectors here
    if(any(temp<2)) cuts[temp<2 ] <- format(cuts[temp<2],justify="left")
    cuts <- paste(cuts, ifelse(temp >=2, ",",
                               ifelse(temp==1, " to the right,", " to the left, ")),
                  sep = '')
  }

  if (is.null(ff$yval2))
    tprint <- x$functions$summary(ff$yval[rows], ff$dev[rows],
                                  ff$wt[rows], ylevel, digits)
  else
    tprint <- x$functions$summary(ff$yval2[rows,,drop=TRUE], ff$dev[rows],
                                  ff$wt[rows], ylevel, digits)


  # Get the nodes
  nodes.id <- match( nodes, id[rows] )

  # Original:   for (ii in 1:length(rows)) {
  for (ii in nodes.id) {
    i <- rows[ii]
    nn <- ff$n[i]
    cat("\nNode number ", id[i], ": ", nn, " observations", sep='')
    if (ff$complexity[i] < cp || is.leaf[i]) cat("\n")
    else cat(",    complexity param=",
             format(signif(ff$complexity[i], digits)), "\n", sep="")

    cat(tprint[ii], "\n")
    if (ff$complexity[i] > cp && !is.leaf[i] ){
      sons <- 2*id[i] + c(0,1)
      sons.n <- ff$n[match(sons, id)]
      cat("  left son=", sons[1], " (", sons.n[1], " obs)",
          " right son=", sons[2], " (", sons.n[2], " obs)", sep='')
      j <- nn - (sons.n[1] + sons.n[2])
      if (j>1) cat(", ", j, " observations remain\n", sep='')
      else if (j==1) cat(", 1 observation remains\n")
      else     cat("\n")
      cat("  Primary splits:\n")
      j <- seq(index[i], length.out=1+ff$ncompete[i])
      if (all(nchar(cuts[j], "w") < 25))
        temp <- format(cuts[j], justify="left")
      else  temp <- cuts[j]
      cat(paste("      ", format(sname[j], justify="left"), " ", temp,
                " improve=", format(signif(x$splits[j,3], digits)),
                ", (", nn - x$splits[j,1], " missing)", sep=''),
          sep="\n")
      if (ff$nsurrogate[i] >0) {
        cat("  Surrogate splits:\n")
        j <- seq(1 +index[i] + ff$ncompete[i], length.out=ff$nsurrogate[i])
        agree <- x$splits[j,3]
        if (all(nchar(cuts[j], "w") < 25))
          temp <- format(cuts[j], justify="left")
        else  temp <- cuts[j]
        if (ncol(x$splits)==5) {
          adj   <- x$splits[j,5]
          cat(paste("      ", format(sname[j], justify="left"), " ",
                    temp,
                    " agree=", format(round(agree, 3)),
                    ", adj=" , format(round(adj, 3)),
                    ", (", x$splits[j,1], " split)", sep=''),
              sep="\n")
        }
        else {                  #an older style rpart object -- no adj value present
          cat(paste("      ", format(sname[j], justify="left"), " ",
                    temp,
                    " agree=", format(round(agree, 3)),
                    ", (", x$splits[j,1], " split)", sep=''),
              sep="\n")
        }
      }
    }
  }
  cat("\n")
  invisible(x)
}

###

# Snip for rpart.plot -------

Snip <- function(...){
  plot(..., snip=TRUE)
}

# Flop  --------------

"Flop" <- function(x, n = 1, m = 1) {

    # flop.rpart(rpart object)
    # Flop branches of rpart object x such that splits are always > (or :), never <
    # B. Compton, 2 March 2003

    #return(x)

    if(n == 1) {
      y <- x
      x <- x$frame
    }

    i <- x[as.numeric(row.names(x)) == n,  ]
    row.names(i) <- as.character(m)

    if(i[1] == "<leaf>")
      i
    else {
      q <- c(2 * n, 2 * n + 1)	# Existing node #s
      r <- c(2 * m, 2 * m + 1)	# New node #s
      print(i$splits)
      if("<" == substring(i$splits[1], 1, 1)) {
        i$splits[1] <- paste(">", substring(i$splits[1], 2),
                             sep = "")
        i$splits[2] <- paste("<", substring(i$splits[2], 2),
                             sep = "")
        q <- rev(q)
      }
      z <- rbind(Flop(x, q[1], r[1]), Flop(x, q[2], r[2])
      )
      z <- rbind(i, z)	# if within recursion, return new node
      if(n != 1)
        z
      else {
        y$frame <- z
        y
      }
    }
  }




# BestTree --------------

BestTree <- function(x) {

  # pick.tree
  # Pick best tree from rpart object using 1 SE rule
  # Returns number of leaves in global leaves
  # B. Compton, 12 May 2004, 21 May 2004

  # example:  BestTree(r.glass)

  i <- x$cptable[,4] <= min(x$cptable[,4] + x$cptable[,5])

  list( leaves = 1 + as.numeric(x$cptable[,2][i][1]),
        cp = as.numeric(x$cptable[,1][i][1] + .000001)   # Add a little for rounding errors
  )
}



# Splits --------------

Splits <- function(x) {

    # original:  splits.rpart
    # B. Compton, 19 Apr 2004
    # Gives results from old z$frame$splits

    # example: Splits(r.glass)

    x <- labels(x, collapse = FALSE)
    x[a <- apply(x, 2, "==", "<leaf>")] <- ""
    a <- (!a) & apply(x, 2, substring, 1, 1) != "<" & apply(x, 2, substring, 1, 1) != ">"
    x[a] <- paste(":", x[a], sep="")
    x[!a] <- paste( apply(x, 2, substr, 1, 1),apply(x, 2, substring, 3, last=100), sep="")[!a]
    colnames(x) <- c("cutleft", "cutright")

    x
}




# LeafRates and LeafHist --------------


#
# # Purity ------------------
# Purity <- function(x, leaves=TRUE) {
#
#   if(x$method != "anova"){
#     # returns the purity of all nodes
#
#     n <- ( ncol(x$frame$yval2) - 2 ) / 2
#     xp <- diag(apply(x$frame$yval2[,-c(1:(n+1))], 1, "[", x$frame$yval))
#
#
#     freq <- data.frame(node=row.names(x$frame), x$frame$yval2[, 2:(1+length(attr(x, "ylevels")))])
#     colnames(freq)[-1] <- attr(x, "ylevels")
#
#     # res <- data.frame(node=rownames(x$frame), purity=xp)
#     names(xp) <- rownames(x$frame)
#
#     if(leaves){
#       xp <- xp[x$frame$var=="<leaf>"]
#       freq <- freq[x$frame$var=="<leaf>",]
#     }
#
#     res <- list(purity=xp, freq=freq)
#
#     class(res) <- c("Purity", class(res))
#     return(res)
#
#   } else {
#     return(NA)
#   }
# }


# Columns of frame include the node, a factor giving the names of the variables used in
# the split at each node (leaf nodes are denoted by the level "<leaf>"), n, the number
# of observations reaching the node, wt, the sum of case weights for observations
# reaching the node, dev, the deviance of the node, yval, the fitted value of the
# response at the node, and splits, a two column matrix of left and right split labels
# for each node. Also included in the frame are complexity, the complexity parameter
# at which this split will collapse, ncompete, the number of competitor splits recorded,
# and nsurrogate, the number of surrogate splits recorded.



LeafRates <- function(x) {

  # original name: rates.rpart
  # Get misclassification rates for rpart object in all end nodes
  # B. Compton, 13 May 2004
  # Rewritten, 16 Nov 2004
  # Modified for root trees, 15 Dec 2004

  if(x$method == "anova"){
    warning("Leafrates can't be returned for regression trees.")
    return(NA)
  }

  xx <- x$frame$yval2[x$frame$var == "<leaf>",]

  if(is.matrix(xx)) {
    z <- matrix(0,dim(xx)[1], 2)
    for(i in 1:dim(xx)[1])
      z[i,1] <- xx[i,xx[i,1]+1]
    z[,2] <- rowSums(xx[,2:((dim(xx)[2]-1)/2+1)]) - z[,1]
  }
  else z <- matrix(xx[c(3,2)],1,2)

  colnames(z) <- c('right', 'wrong')
  rownames(z) <- rownames(x$frame[x$frame$var == "<leaf>",])


  structure(
    list(node=rownames(z),
         freq=z,
         p.row=prop.table(z, 1),
         mfreq=apply(z, 1, sum),
         mperc=prop.table(apply(z, 1, sum))
         ),
    class=c("LeafRates","list"))

}



print.LeafRates <- function(x, ...){

  z <- rbind(c("","freq","","perc","","total",""),
               c("node","right","wrong","right","wrong","abs","perc"),
               "",
               cbind(x$node,
                     Format(x$freq, fmt="abs"),
                     Format(x$p.row, fmt="%", digits=1),
                     Format(x$mfreq, fmt="abs"),
                     Format(x$mperc, fmt="%")))

  write.table(format(z, justify="left"), sep = " ",
              row.names=FALSE, col.names=FALSE, quote=FALSE)

}





plot.LeafRates <- function(x, col=NULL, which=c("rel","abs"),
                           layout=NULL, ylim=NULL, ...){



  # m <- x$frame$yval2[,-1]
  # m <- m[,1:((ncol(m)-1)/2)]
  # rownames(m) <- rownames(x$frame)
  # m.leaves <- m[x$frame$var == "<leaf>",]
  #
  # m.leaves <- x
  # if(type=="rel") m.leaves <- prop.table(m.leaves, 1)
  # d.leaves <- data.frame(frq=as.vector(t(m.leaves)),
  #                        leaf=factor(rep(rownames(m.leaves), each=length(attr(x, "ylevels"))),
  #                                    levels=rep(rownames(m.leaves))),
  #                        class=rep(attr(x, "ylevels"), times=nrow(m.leaves))
  # )

  if(is.null(col))
    col <- c(DescTools::hblue, DescTools::hred)

  which <- match.arg(which)

  if(which=="rel"){
    if(is.null(ylim))
      ylim = c(0, 1.1)
    ylab = "relative freq."

  } else {
    ylab = "frequency"

  }

  d.leaves <- data.frame(frq=c(x$freq[,1], x$freq[,2]),
                           leaf=rep(row.names(x$freq), times=2),
                           # class=rep(colnames(x$freq), each=nrow(x$freq)), stringsAsFactors = FALSE)
                         class=factor(rep(c("T","F"), each=nrow(x$freq)), levels=c("T","F")), stringsAsFactors = FALSE)

  idx <- d.leaves$leaf[1:nrow(x$freq)]

  if(which=="rel")
    d.leaves$frq <- c(x$p.row[,1], x$p.row[,2])

  if(is.null(layout))
    layout <- c(nrow(x$freq), 1)

  # what a zangengeburt.... node order as in the original result
  id <- order(as.numeric(Sort(cbind(nod=rownames(x$freq),
                                    nr=1:nrow(x$freq)))[,2]))

  barchart( frq ~ class|leaf,  data=d.leaves,
            stack=FALSE, horiz=FALSE, layout=layout, ylab=ylab, xlab="classes",
            ylim=ylim, index.cond=list(id),
            panel = function(y, x, ...){
              cols <- rep(col, length(y))
              cols[which.max(y)] <- col[1]
              panel.barchart(x, y,..., col=cols, origin=0, box.width=0.85)
            },

            strip=.myStripStyle
  )

}



# Complexity parameter
CP <- function(x, ...){
  structure(list(cp=x$cptable,
                 mincp=x$cptable[which.min(x$cptable[,"xerror"]), "CP"],
                 x=x), class=c("CP"))
}



print.CP <- function(x, digits = getOption("digits") - 2L, ...){
  rpart::printcp(x$x, digits=digits, ...)
  cat("\n")
}

plot.CP <- function(x, minline = TRUE, lty = 3,
                    col = 1, upper = c("size", "splits", "none"), ...){
  rpart::plotcp(x$x, minline, lty, col, upper=upper, ...)
}





# plot.rpart ------------------

# we override the native rpart plot and use rpart.plot here

plot.rpart <- function (x = stop("no 'x' arg"), type = 2, extra = "auto",
                        under = FALSE, fallen.leaves = TRUE, digits = 2, varlen = 0,
                        faclen = 0, roundint = TRUE, cex = NULL, tweak = 1, clip.facs = FALSE,
                        clip.right.labs = TRUE, snip = FALSE, box.palette = "auto",
                        shadow.col = 0, node.labels=TRUE, ...) {

  if(identical(box.palette, "auto"))
    box.palette <- as.list(ColToOpaque(SetAlpha(Pal("Helsana"))))

  b <- rpart.plot(x =x, type = type, extra = extra,
                  under = under, fallen.leaves = fallen.leaves, digits = digits, varlen = varlen,
                  faclen = faclen, roundint = roundint, cex = cex, tweak = tweak,
                  clip.facs = clip.facs,
                  clip.right.labs = clip.right.labs, snip = snip, box.palette =  box.palette ,
                  shadow.col = shadow.col, ...)

  if(node.labels){
    par(xpd=TRUE)
    BoxedText(x = b$boxes$x1, y=b$boxes$y2, labels = rownames(x$frame),
              border=NA, txt.col = "steelblue", font=2, col=SetAlpha("white", 0.7), cex=0.8)
  }

  invisible(b)

}



.predict.rpart <- function (object, newdata, type = c("vector", "prob", "class",
                                   "matrix", "where", "leaf"), na.action = na.pass, ...){

  type <- match.arg(type)

  if(type %in% c("where","leaf"))
    .predict.leaves(object, newdata, type = type)
  else {
    class(object) <- class(object)[class(object) != "FitMod"]
    predict(object, newdata, type, na.action, ...)
  }


}


.predict.leaves <- function (rp, newdata, type = "where") {

  if (type == "where") {
    rp$frame$yval <- 1:nrow(rp$frame)
    should.be.leaves <- which(rp$frame[, 1] == "<leaf>")
  }
  else if (type == "leaf") {
    rp$frame$yval <- rownames(rp$frame)
    should.be.leaves <- rownames(rp$frame)[rp$frame[, 1] ==
                                             "<leaf>"]
  }
  else stop("Type must be 'where' or 'leaf'")
  leaves <- predict(rp, newdata = newdata, type = "vector")
  should.be.leaves <- which(rp$frame[, 1] == "<leaf>")
  bad.leaves <- leaves[!is.element(leaves, should.be.leaves)]
  if (length(bad.leaves) == 0)
    return(leaves)
  u.bad.leaves <- unique(bad.leaves)
  u.bad.nodes <- row.names(rp$frame)[u.bad.leaves]
  all.nodes <- row.names(rp$frame)[rp$frame[, 1] == "<leaf>"]
  is.descendant <- function(all.leaves, node) {
    if (length(all.leaves) == 0)
      return(logical(0))
    all.leaves <- as.numeric(all.leaves)
    node <- as.numeric(node)
    if (missing(node))
      return(NA)
    result <- rep(FALSE, length(all.leaves))
    for (i in 1:length(all.leaves)) {
      LEAF <- all.leaves[i]
      while (LEAF > node) {
        LEAF <- trunc(LEAF/2)
        if (LEAF == node)
          result[i] <- TRUE
        break
      }
    }
    return(result)
  }
  where.tbl <- table(rp$where)
  names(where.tbl) <- row.names(rp$frame)[as.numeric(names(where.tbl))]
  for (u in 1:length(u.bad.nodes)) {
    desc.vec <- is.descendant(all.nodes, u.bad.nodes[u])
    me <- where.tbl[all.nodes][desc.vec]
    winner <- names(me)[me == max(me)][1]
    leaves[leaves == u.bad.leaves[u]] <- which(row.names(rp$frame) ==
                                                 winner)
  }
  return(leaves)
}




Node <- function (x, node=NULL, type=c("all", "split", "leaf"), digits=3) {

  if (!inherits(x, "rpart"))
    stop("Not a legitimate \"rpart\" object")

  ff <- x$frame
  ylevel <- attr(x, "ylevels")
  id <- as.integer(row.names(ff))

  parent.id <- ifelse(id == 1L, 1L, id%/%2L)

  rows <- seq_along(id)
  is.leaf <- ff$var == "<leaf>"
  index <- cumsum(c(1L, ff$ncompete + ff$nsurrogate + !is.leaf))

  if (!all(is.leaf)) {
    sname <- rownames(x$splits)
    cuts <- character(nrow(x$splits))
    temp <- x$splits[, 2L]
    for (i in seq_along(cuts)) {
      cuts[i] <- if (temp[i] == -1L)
        paste("<", format(signif(x$splits[i, 4L], digits)))
      else if (temp[i] == 1L)
        paste("<", format(signif(x$splits[i, 4L], digits)))
      else paste("splits as ", paste(c("L", "-", "R")[x$csplit[x$splits[i, 4L], 1:temp[i]]],
                                     collapse = "", sep = ""), collapse = "")
    }
    if (any(temp < 2L))
      cuts[temp < 2L] <- format(cuts[temp < 2L], justify = "left")

    cuts <- paste0(cuts, ifelse(temp >= 2L, ",", ifelse(temp == 1L, " to the right,", " to the left, ")))
  }

  tmp <- if (is.null(ff$yval2))
    ff$yval[rows]
  else
    ff$yval2[rows, , drop = FALSE]
  tmp <- unname(tmp)

  tprint <- x$functions$summary(tmp, ff$dev[rows], ff$wt[rows], ylevel, digits)
  nclass <- (ncol(tmp) - 2L)/2L

  nlst <- list()

  if(is.null(node))
    nid <- seq_along(rows)
  else
    nid <- seq_along(rows)[which(!is.na(rows[match(row.names(ff), node)]))]

  for (ii in nid) {

    i <- rows[ii]
    nlbl <- as.character(id[i])

    nlst[[nlbl]] <- SetNames(list(), names=nlbl)
    nn <- ff$n[i]
    nlst[[nlbl]]$id <- id[i]
    nlst[[nlbl]]$vname <- as.character(ff$var[i])
    nlst[[nlbl]]$isleaf <- is.leaf[i]
    nlst[[nlbl]]$nobs <- nn

    nlst[[nlbl]]$group <- ylevel[tmp[i, 1L]]
    nlst[[nlbl]]$ycount <- tmp[i, 1L + (1L:nclass)]
    nlst[[nlbl]]$yprob <- tmp[i, 1L + nclass + 1L:nclass]
    nlst[[nlbl]]$nodeprob <- tmp[i, 2L * nclass + 2L]

    nlst[[nlbl]]$complexity <- ff$complexity[i]

    nlst[[nlbl]]$tprint <- tprint[ii]

    if (!is.leaf[i]) {
      sons <- 2L * id[i] + c(0L, 1L)
      sons_n <- ff$n[match(sons, id)]

      nlst[[nlbl]]$sons <- SetNames(sons, names=c("left", "right"))
      nlst[[nlbl]]$sons_n <- SetNames(sons_n, names=c("left", "right"))


      ## what's that??
      # j <- nn - (sons.n[1L] + sons.n[2L])
      # if (j > 1L)
      #   cat(", ", j, " observations remain\n",
      #       sep = "")
      # else if (j == 1L)
      #   cat(", 1 observation remains\n")
      # else cat("\n")


      # primary splits
      j <- seq(index[i], length.out = 1L + ff$ncompete[i])
      temp <- if (all(nchar(cuts[j], "w") < 25L))
        format(cuts[j], justify = "left")
      else
        cuts[j]

      nlst[[nlbl]]$primarysplits <- data.frame(split=sname[j], direction=temp, improve=x$splits[j, 3L], missing= nn - x$splits[j, 1L])

      # surrogate splits
      if (ff$nsurrogate[i] > 0L) {
        j <- seq(1L + index[i] + ff$ncompete[i], length.out = ff$nsurrogate[i])
        agree <- x$splits[j, 3L]
        temp <- if (all(nchar(cuts[j], "w") < 25L))
          format(cuts[j], justify = "left")
        else cuts[j]
        adj <- x$splits[j, 5L]

        nlst[[nlbl]]$surrogatesplits <- data.frame(split=sname[j], direction=temp,
                                                   agree=agree, adj=adj, split=x$splits[j, 1L])

      }
    }
  }

  type <- match.arg(type)
  if(type=="leaf")
    nlst <- nlst[sapply(nlst, "[[", "isleaf")]

  if(type=="split")
    nlst <- nlst[!sapply(nlst, "[[", "isleaf")]


  return(structure(nlst, class="node"))

}


print.node <- function(x, digits=3, ...){

  if(length(x)==0)
    cat("list()\n")

  else {

    for(i in seq_along(x)) {
      cat("\nNode number ", x[[i]]$id, ": ", x[[i]]$nobs, " observations", sep = "")

      if (x[[i]]$isleaf)
        cat("\n")
      else
        cat(",    complexity param=", format(signif(x[[i]]$complexity, digits)), "\n", sep = "")

      cat(x[[i]]$tprint, "\n")

      if(!x[[i]]$isleaf){
        cat("  left son=", x[[i]]$sons[1L], " (", x[[i]]$sons_n[1L],
            " obs)", " right son=", x[[i]]$sons[2L],
            " (", x[[i]]$sons_n[2L], " obs)", sep = "")


        cat("\n  Primary splits:\n")

        for(lx in split(x[[i]][["primarysplits"]],
                        f = 1:nrow(x[[i]][["primarysplits"]]))) {

          cat(paste("      ", format(lx$split, justify = "left"),
                    " ", lx$direction, " improve=", format(signif(lx$improve, digits)),
                    ", (",
                    lx$missing, " missing)", sep = ""), sep = "\n")

        }

        if(!is.null(x[[i]]$surrogatesplits)){
          cat("  Surrogate splits:\n")
          for(lx in split(x[[i]][["surrogatesplits"]],
                          f = 1:nrow(x[[i]][["surrogatesplits"]]))) {

            cat(paste("      ", format(lx$split, justify = "left"),
                      " ", lx$direction, " agree=", format(round(lx$agree, 3L)),
                      ", adj=", format(round(lx$adj, 3L)),
                      ", (", lx$split.1, " split)",
                      sep = ""), sep = "\n")

          }
        }
      }

      cat("\n")

    }

  }
}



#' #' @importFrom flipU OutcomeName
#' #' @importFrom stats formula
#' rPartToTreeFrame <- function(obj)
#' {
#'   frame <- obj$frame
#'   n.nodes <- nrow(frame)
#'   splits <- matrix("", n.nodes, 2)
#'   nms <- colnames(obj$splits)
#'   index.i <- nms == "index"
#'   ncat.i <- nms == "ncat"
#'   c <- 1
#'   for (i in 1:n.nodes)
#'   {
#'     var.name <- as.character(frame$var[i])
#'     if (var.name != "<leaf>")
#'     {
#'       if (is.factor(obj$model[[var.name]]))
#'       {
#'         n.levels <- length(attr(obj, "xlevels")[[var.name]])
#'         extended.letters <- extendLetters(n.levels)
#'         lttrs <- extended.letters[1:n.levels]
#'         levels.split <- obj$csplit[obj$splits[c, index.i], ]
#'         splits[i, 1] <- paste0(":", paste(lttrs[levels.split == 1], collapse = ""))
#'         splits[i, 2] <- paste0(":", paste(lttrs[levels.split == 3], collapse = ""))
#'       }
#'       else
#'       {
#'         break.val <- obj$splits[c, index.i]
#'         if (obj$splits[c, ncat.i] < 0)
#'         {
#'           splits[i, 1] <- paste0("<", break.val)
#'           splits[i, 2] <- paste0(">=", break.val)
#'         }
#'         else
#'         {
#'           splits[i, 1] <- paste0(">=", break.val)
#'           splits[i, 2] <- paste0("<", break.val)
#'         }
#'       }
#'       c <- c + frame$ncompete[i] + frame$nsurrogate[i] + 1
#'     }
#'   }
#'   frame$splits <- splits
#'
#'   outcome.var <- obj$model[[OutcomeName(formula(obj$terms))]]
#'
#'   if (is.factor(outcome.var))
#'     frame$yval <- levels(outcome.var)[frame$yval]
#'
#'   if (!is.null(frame$yval2))
#'   {
#'     fitted.levels <- levels(droplevels(outcome.var[obj$subset]))
#'     all.levels <- levels(outcome.var)
#'     n.levels <- length(fitted.levels)
#'     yprob <- frame$yval2[, (n.levels + 2):(2 * n.levels + 1)]
#'     empty.levels <- all.levels[!all.levels %in% fitted.levels]
#'     yprob <- cbind(yprob, matrix(0, nrow = nrow(yprob), ncol = length(empty.levels)))
#'     colnames(yprob) <- c(fitted.levels, empty.levels)
#'     frame$yprob <- yprob
#'   }
#'   frame
#' }


