.First.lib <- function(library, pkg)
{
    library.dynam("sm", pkg, library)
    assign(".sm.home", file.path(library, pkg),
           pos=match("package:sm", search()))
  ".sm.Options" <- assign(".sm.Options",
      list(
        hmult=1, h.weights=NA, describe=TRUE, diff.ord=2,
        add=FALSE,  band=TRUE,  props = c(75,50,25),
        nbins=NA, positive=FALSE,  delta=NA, display=NA,
        xlab = NA, ylab = NA, zlab = NA, xlim = NA, ylim = NA, yht = NA,
        panel=FALSE, ngrid = NA, eval.points = NA, rugplot=TRUE,
        lty=1, col=1,  pch=1 , theta=pi/4, phi=pi/4,
        poly.index=1, test=TRUE, eye.mult=c(-6, -8, 5), hull=TRUE),
        pos=match("package:sm", search()) )
  cat("Library `sm', version 2; Copyright (C) 1997, 2000 A.W.Bowman & A.Azzalini\n")
  cat("type help(sm) for summary information\n")
    invisible()
}

isMatrix <- function(x) length(dim(x)) == 2

isInteger <- function(x) all(x == round(x))

sm.script <- function(name, path)
{
    if (missing(path)) path <- file.path(.sm.home, "scripts")
    else path <- as.character(substitute(path))
    if (missing(name)) {
      file.show(file.path(path, "index.doc"))
    } else {
        name <- as.character(substitute(name))
        if(length(name) == 3 && name[1] == "<-")
            name <- paste(name[2:3], collapse="_")
        source(file.path(path, paste(name, ".q", sep = "")))
    }
    invisible()
}
