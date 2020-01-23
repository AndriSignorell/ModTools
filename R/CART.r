# ****************************************************************************
#
# Projekt:  RegClassTools.r (will be sourced there)
#
# Zweck:    some CART stuff, routines for handling rparts
#
# Autor:    partly by B. Compton, compiled and extended by Andri Signorell
# Version:	0.1 (in development)
#
# Depends:
#
# Datum:
#           04.12.2013 	created
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



# Rules ------------------
Rules <- function(x, node=NULL, leafonly = FALSE) {

  if (!inherits(x, "rpart")) stop("Not a legitimate rpart tree")
  #
  # Get some information.
  #

  node <- as.character(node)

  if(is.null(node))
    node <- rownames(x$frame)

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




LeafRates <- function(x) {

  # original name: rates.rpart
  # Get misclassification rates for rpart object in all end nodes
  # B. Compton, 13 May 2004
  # Rewritten, 16 Nov 2004
  # Modified for root trees, 15 Dec 2004

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
    list(freq=z, p.row=prop.table(z, 1), mfreq=apply(z, 1, sum), mperc=prop.table(apply(z, 1, sum))),
    class=c("LeafRates","list"))

}


print.LeafRates <- function(x, ...){

  tlst <- list(
    freq = Format(x$freq, fmt=Fmt("abs")),
    perc = Format(x$p.row, fmt=Fmt("per")),
    total = cbind(Format(x$mfreq, fmt=Fmt("abs")), Format(x$mperc, fmt=Fmt("per")))
  )

  ma <- do.call("Abind", c(tlst, along = 3, use.dnns = TRUE))
  ft <- ftable(ma, col.vars=c(3,2))

  print(ft, justify="right", ...)

}





plot.LeafRates <- function(x, col=c(hblue, hred), type=c("rel","abs"),
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

  type <- match.arg(type)

  if(type=="rel"){
    if(is.null(ylim))
      ylim = c(0, 1.1)
    ylab = "relative freq."

  } else {
    ylab = "frequency"

  }

  d.leaves <- data.frame(frq=c(x$freq[,1], x$freq[,2]),
                           leaf=rep(row.names(x$freq), times=2),
                           class=rep(colnames(x$freq), each=nrow(x$freq)), stringsAsFactors = FALSE)

  idx <- d.leaves$leaf[1:nrow(x$freq)]

  if(type=="rel")
    d.leaves$frq <- c(x$p.row[,1], x$p.row[,2])

  if(is.null(layout))
    layout <- c(nrow(x$freq), 1)

  # what a zangengeburt.... node order as in the original result
  id <- order(as.numeric(Sort(cbind(nod=rownames(x$freq), nr=1:nrow(x$freq)))[,2]))

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



.myStripStyle <- function(which.panel, factor.levels, ...) {
  panel.rect(0, -0.5, 1, 1,
             col = "grey",
             border = 1)
  panel.text(x = 0.5, y = 0.25,
             font=2,
             lab = gettextf("node %s", factor.levels[which.panel]),
             col = "black")
}





# Purity ------------------
Purity <- function(x, leaves=TRUE) {

  if(x$method != "anova"){
    # returns the purity of all nodes

    n <- ( ncol(x$frame$yval2) - 2 ) / 2
    xp <- diag(apply(x$frame$yval2[,-c(1:(n+1))], 1, "[", x$frame$yval))


    freq <- data.frame(node=row.names(x$frame), x$frame$yval2[, 2:(1+length(attr(x, "ylevels")))])
    colnames(freq)[-1] <- attr(x, "ylevels")

    # res <- data.frame(node=rownames(x$frame), purity=xp)
    names(xp) <- rownames(x$frame)

    if(leaves){
      xp <- xp[x$frame$var=="<leaf>"]
      freq <- freq[x$frame$var=="<leaf>",]
    }

    res <- list(purity=xp, freq=freq)

    class(res) <- c("Purity", class(res))
    return(res)

  } else {
    return(NA)
  }
}




plot.Purity <- function(x, col=c(hblue, hred), type=c("abs","rel"),
                        layout=NULL, ylim=NULL, ...){

  # myStripStyle <- function(which.panel, factor.levels, ...) {
  #   panel.rect(0, -0.5, 1, 1,
  #              col = "grey",
  #              border = 1)
  #   panel.text(x = 0.5, y = 0.25,
  #              font=2,
  #              lab = gettextf("%s", factor.levels[which.panel]),
  #              col = "black")
  # }


  ylab <- ""

  if(is.null(layout))
    layout <- c(nrow(x$freq), 1)


  d.nodes <- reshape(x$freq, idvar="node", varying=c("neg","pos"), v.names="n", timevar="ind",
          direction="long")

  if(is.null(ylim))
    ylim <- c(0, max(pretty(d.nodes$n*1.1)))

  # what a zangengeburt.... node order as in the original result
  id <- order(as.numeric(Sort(cbind(nod=as.character(x$freq$node), nr=1:nrow(x$freq)))[,2]))

  barchart( n ~ ind | node,  data=d.nodes,
            stack=FALSE, horiz=FALSE, layout=layout, ylab=ylab, xlab="classes",
            ylim=ylim, index.cond=list(id),
            panel = function(y, x, ...){
              cols <- rep(col, length(y))
              # cols[which.max(y)] <- col[1]
              panel.barchart(x, y,..., col=cols, origin=0, box.width=0.85)
            },

            strip=.myStripStyle
  )

}





# plot.rpart ------------------
# we override the native rpart plot
plot.rpart <- rpart.plot


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





PlotTree <- function(x, type=c("vert", "horiz", "min"), cols=NULL, shade=NULL,
                     args.legend=NULL, main="", fallen.leaves = TRUE, ...){

  # example:
  # PlotTree(r.rpart, cols=cols, uniform=TRUE, args.legend=list(x="topleft", inset=0.05), main="Type ~ .")

  .nodeFullVert <- function(x, labs, digits, varlen) {
    # function for rpart.plot to describe the leaf nodes with abs. and rel freq
    m <- x$frame$yval2[,-1]
    m <- m[,1:((ncol(m)-1)/2)]
    pm <- gsub(pattern=" 0", replacement=" ",
               apply(apply(apply(m, 1, prop.table),1, Format, digits=3), 1, paste, collapse=" ")
    )

    nclass <- attr(x, "ylevels")[x$frame$yval]
    nnode <- rownames(x$frame)
    txt <- character(nrow(m))
    for(i in 1: nrow(m)){
      ns <- data.frame(class=attr(x, "ylevels"),
                       obs=m[i,], perc=Format(t(apply(m, 1, prop.table))[i,] *100, 1), stringsAsFactors=FALSE)
      ns$class[which.max(ns$perc)] <- paste(">", ns$class[which.max(ns$perc)])
      txt[i] <- paste(capture.output(print(ns, row.names=FALSE, digits=digits)), collapse="\n")
      txt[i] <- paste(gettextf("  node %s - class %s\n     n = %s\n-------------------\n",
                               nnode[i], nclass[i], x$frame$n[i]), txt[i])
    }
    return(txt)
  }

  .nodeFull <- function(x, labs, digits, varlen) {
    # function for rpart.plot to describe the leaf nodes with abs. and rel freq
    m <- x$frame$yval2[,-1]
    m <- m[,1:((ncol(m)-1)/2)]
    pm <- gsub(pattern=" 0", replacement=" ",
               apply(apply(apply(m, 1, prop.table),1, Format, digits=3), 1, paste, collapse=" ")
    )

    paste(gettextf("node %s, class %s\n", rownames(x$frame), labs), apply(m, 1, paste, collapse=" "),
          "\n", pm)
  }

  .nodeNone <- function(x, labs, digits, varlen){rep("  ", nrow(x$frame))}
  .splitNone <- function(x, labs, digits, varlen, faclen){rep("  ", nrow(x$frame))}

  xpd <- par(xpd=TRUE); on.exit(par(xpd))

  if(is.null(cols)){
    cols <- "white"
    if(is.null(args.legend)) args.legend <- NA  # no default legend if there are no colors
  }
  nodecols <- cols[x$frame$yval2[,1]]

  if(is.null(shade)) shade <- TRUE
  if(shade & x$method=="class"){
    # we don't have purity for regression trees
    nodecols <- SetAlpha(nodecols, Purity(x)$purity)
  }

  switch( match.arg(type) ,
          "vert" = {
            res <- plot(x, node.fun=.nodeFullVert, box.col=nodecols, fallen.leaves=fallen.leaves,
                        split.box.col="lightgray",   # lightgray split boxes (default is white)
                        split.border.col="darkgray", # darkgray border on split boxes
                        split.round=.4,
                        leaf.round=.4, main = main,
                        ...)
          },
          "horiz" = {
            res <- plot(x, node.fun=.nodeFull, box.col=nodecols, fallen.leaves=fallen.leaves,
                        split.box.col="lightgray",   # lightgray split boxes (default is white)
                        split.border.col="darkgray", # darkgray border on split boxes
                        split.round=.4,
                        leaf.round=.4, main = main,
                        ...)
          },
          "min" =  {
            res <- plot(x, node.fun=.nodeNone, box.col=nodecols, round=0
                        , split.shadow.col="grey", split.box.col=nodecols, split.border.col="black"
                        , split.fun=.splitNone, yesno=FALSE
                        , boxes.include.gap = TRUE, main=main, nspace=0, minbranch=.01, ...)
          }
  )

  args.legend1 <- list( x="topright", inset=0.02, legend=attr(x, "ylevels")
                        , fill=col, bg="white", cex=0.8 )
  args.legend1[["fill"]] <- cols
  args.legend1[["title"]] <- " classes"

  if ( !is.null(args.legend) ) { args.legend1[names(args.legend)] <- args.legend }
  add.legend <- TRUE
  if(!is.null(args.legend)) if(all(is.na(args.legend))) {add.legend <- FALSE}

  if(add.legend) do.call("legend", args.legend1)

  invisible(res)

}

# PlotTree(r.pima, uniform=FALSE, type="horiz", fallen.leaves=TRUE)

# PlotTree(r.rp, uniform=FALSE, type="horiz", fallen.leaves=TRUE)

#
# PlotTree examples:
#
# windows(w=20, h=15)
# PlotTree(r.glass, type="fullv", cols=PalTibco()[2:8], fallen.leaves=FALSE,
#          args.legend=list(x="topleft"), main="Glass")
#
# windows(w=6, h=5)
# PlotTree(r.glass, type="min", cols=PalTibco()[2:8], main="Glass",
#          fallen.leaves=FALSE, args.legend=list(x="topleft"))

# PlotTree(r.rpart, type="fullv", cols=cols, uniform=TRUE, args.legend=list(x="topleft", inset=0.05), main="Type ~ .")





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



Node <- function(x, node, digits = getOption("digits")){

  node <- as.character(node)

  cp  <-  0

  ff <- x$frame
  ylevel <- attr(x, "ylevels")
  id <- as.integer(row.names(ff))
  parent.id <- ifelse(id == 1L, 1L, id%/%2L)
  parent.cp <- ff$complexity[match(parent.id, id)]
  rows <- seq_along(id)[parent.cp > cp]
  rows <- if (length(rows))
    rows[order(id[rows])]
  else 1L
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
      else paste("splits as ", paste(c("L", "-", "R")[x$csplit[x$splits[i,
                                                                        4L], 1:temp[i]]], collapse = "", sep = ""), collapse = "")
    }
    if (any(temp < 2L))
      cuts[temp < 2L] <- format(cuts[temp < 2L], justify = "left")
    cuts <- paste0(cuts, ifelse(temp >= 2L, ",", ifelse(temp ==
                                                          1L, " to the right,", " to the left, ")))
  }
  tmp <- if (is.null(ff$yval2))
    ff$yval[rows]
  else ff$yval2[rows, , drop = FALSE]
  tprint <- x$functions$summary(tmp, ff$dev[rows], ff$wt[rows],
                                ylevel, digits)
#  for (ii in seq_along(rows)) {
  for (ii in which(row.names(ff)[rows] %in% node)) {
    i <- rows[ii]
    nn <- ff$n[i]
    cat("\nNode number ", id[i], ": ", nn, " observations",
        sep = "")
    if (ff$complexity[i] < cp || is.leaf[i])
      cat("\n")
    else cat(",    complexity param=", format(signif(ff$complexity[i],
                                                     digits)), "\n", sep = "")
    cat(tprint[ii], "\n")
    if (ff$complexity[i] > cp && !is.leaf[i]) {
      sons <- 2L * id[i] + c(0L, 1L)
      sons.n <- ff$n[match(sons, id)]
      cat("  left son=", sons[1L], " (", sons.n[1L], " obs)",
          " right son=", sons[2L], " (", sons.n[2L], " obs)",
          sep = "")
      j <- nn - (sons.n[1L] + sons.n[2L])
      if (j > 1L)
        cat(", ", j, " observations remain\n", sep = "")
      else if (j == 1L)
        cat(", 1 observation remains\n")
      else cat("\n")
      cat("  Primary splits:\n")
      j <- seq(index[i], length.out = 1L + ff$ncompete[i])
      temp <- if (all(nchar(cuts[j], "w") < 25L))
        format(cuts[j], justify = "left")
      else cuts[j]
      cat(paste("      ", format(sname[j], justify = "left"),
                " ", temp, " improve=", format(signif(x$splits[j,
                                                               3L], digits)), ", (", nn - x$splits[j, 1L],
                " missing)", sep = ""), sep = "\n")
      if (ff$nsurrogate[i] > 0L) {
        cat("  Surrogate splits:\n")
        j <- seq(1L + index[i] + ff$ncompete[i], length.out = ff$nsurrogate[i])
        agree <- x$splits[j, 3L]
        temp <- if (all(nchar(cuts[j], "w") < 25L))
          format(cuts[j], justify = "left")
        else cuts[j]
        adj <- x$splits[j, 5L]
        cat(paste("      ", format(sname[j], justify = "left"),
                  " ", temp, " agree=", format(round(agree, 3L)),
                  ", adj=", format(round(adj, 3L)), ", (", x$splits[j,
                                                                    1L], " split)", sep = ""), sep = "\n")
      }
    }
  }


  cat("\n")
  invisible(x)
}


###



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




# usual --------------

# These are B. Comptons defaults
# "usual" <-
#   structure(list(minsplit = 8, minbucket = 5, cp = 0, maxcompete = 100,
#                  maxsurrogate = 100, usesurrogate = 2, surrogatestyle = 0,
#                  xval = 10, maxdepth = 30), .Names = c("minsplit", "minbucket",
#                                                        "cp", "maxcompete", "maxsurrogate", "usesurrogate", "surrogatestyle",
#                                                        "xval", "maxdepth"))
#

# loss ??? what is this used for?  --------------

# "loss" <-
#   function(cost)
#   {
#     # Create 2x2 loss matrix from cost
#
#     if (length(cost) == 1)
#       cost <- c(cost,1)
#     list(loss = t(matrix(c(0,cost,0),2,2)))
#   }
#


# kappatau
#
# "kappatau" <-
#   function(c, p)
#   {
#     # kappatau
#     # Compute Kappa or Tau statistic for confusion matrix.  Returns as global kt.
#     # Reference: McGarigal, K, S. Cushman, and S. Stafford. 2000. Multivariate statistics
#     #    for wildlife and ecology research. Springer-Verlag, New York. Pages 165-166.
#     # B. Compton, 15 Nov 2004, 12 Aug 2005, 29 Nov 2005, 8 May 2006
#
#     b <- 0 != rowSums(c)+colSums(c)		# remove any unused classes
#     c <- c[b,b]
#     p <- p[b]
#     if (all(p == (rowSums(c)/sum(c)))) {
#       # Priors are proportional to group sizes (default), so do kappa
#       e <- sum(rowSums(c)*colSums(c))/sum(c)
#       kt <<- (sum(diag(c)) - e) / (sum(c) - e)
#
#       cat ('Kappa = ',round(kt,3), '\n')
#     }
#     else {
#       # Priors specified, so do tau
#       e <- sum(p * rowSums(c))
#       kt <<- (sum(diag(c)) - e) / (sum(c) - e)
#
#       cat ('Tau = ',round(kt,3), '\n')
#     }
#   }

# Improve: don't know what this does... ----------------

# Improve <- function(n, l, r, t, c, p, y) {
#
#   .SSD <- function(x) {
#     # Give sum of squared deviations of x
#     sum((x - mean(x))^2)
#   }
#
#   .gini <- function(x, c = NULL, p = NULL, y) {
#
#     # Gini index of x
#     # Use costs c and priors p, and entire response variable y
#     # B. Compton, 29 Nov 2005
#     #
#     # Current version (14 Dec 2005) implements priors but not costs.
#     # It matches Kevin's Gini/priors results, but still doesn't
#     # give me same improvement as rpart for primary splitters when
#     # priors are set (it matches when priors come from data).
#     # Numbers are close (within 15%), but with no clear pattern of
#     # disagreement.  Either Kevin & I have both misinterperted how
#     # the Gini index incorporates priors, or the weightings for
#     # the left/right splits are wrong.  Grr.
#
#
#     g <- length(u <- unique(y)) # number of groups (get levels from response variable)
#     h <- table(y) # total tree count in each group
#     n <- rowSums(outer(u,x,FUN='==')) # number in each group
#     s <- sum(n) # n across groups
#
#     if(is.null(c))    c <- 1*outer(1:g,1:g,FUN='!=') # default costs = all 1
#     if(is.null(p))    p <- h/sum(h) # default priors = from data
#
#     r <- p*n/h                    # Kevin's version
#     1 - sum((r/sum(r))^2)
#
#   }
#
#
#   # Calculate delta-I for values in node, left child, and right child
#   # Use deviance for regression trees (t = 'anova') and gini for classification trees
#   # c = costs, p = priors, y = values of response variable
#
#   catn <- function(...){cat(paste(...,'\n'))}
#
#   n <- n[!is.na(n)]
#   l <- l[!is.na(l)]
#   r <- r[!is.na(r)]
#
#   if(t == 'anova') {  	# If regression tree, use deviance
#     z <- .SSD(n) - .SSD(l) - .SSD(r)
#     if(is.na(z)) catn(n,l,r)
#   }
#   else {					# For classification trees, use Gini index
#
#     catn('length(l)=',length(l),'length(r)=',length(r),'length(n)=',length(n))
#
#     q <- c(length(l),length(r)) / length(n)
#
#     a <- .gini(n, c, p, y)
#     b <- .gini(l, c, p, y)
#     c <- .gini(r, c, p, y)
#
#     z <- a - sum(q * c(b,c))
#     catn('a=',a,'b=',b,'c=',c,'q=',q,'improvement=',z)
#
#   }
#
#   z
# }
