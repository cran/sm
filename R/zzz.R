".sm.Options" <-
    list(hmult = 1, h.weights = NA, describe = TRUE, diff.ord = 2,
         add = FALSE, band = TRUE, props = c(75,50,25), nbins = NA,
         positive = FALSE, delta = NA, display = NA, xlab = NA, 
         ylab = NA, zlab = NA, xlim = NA, ylim = NA, zlim=NA, yht = NA,
         panel = FALSE, ngrid = NA, eval.points = NA, rugplot = TRUE, 
         lty = 1, col = 1, pch = 1 , theta = -30, phi = 40,
         poly.index = 1, test = TRUE, hull = TRUE, verbose=1, df
         = NA, method = NA, structure.2d = "scaled", nboot=100, df=4,
         show.script=TRUE, eval.grid=TRUE)

.onAttach <- function(library, pkg)
{
    ## we can't do this in .onLoad
    unlockBinding(".sm.Options", asNamespace("sm"))
    cat("Library `sm', version 2.1; ",
        "Copyright (C) 1997, 2000, 2005 A.W.Bowman & A.Azzalini\n")
    cat("type help(sm) for summary information\n")
    invisible()
}


isMatrix <- function(x) length(dim(x)) == 2

isInteger <- function(x) all(x == round(x))

sm.script <- function(name, path)
{
    if (missing(path)) path <- system.file("scripts", package = "sm")
    else path <- as.character(substitute(path))
    if (missing(name)) {
        file.show(file.path(path, "index.doc"))
    } else {
        name <- as.character(substitute(name))
        if(length(name) == 3 && name[1] == "<-")
            name <- paste(name[2:3], collapse="_")
        file <- file.path(path, paste(name, ".q", sep = ""))
        if(sm.options()$show.script)
          file.show(file, title=name, header=paste('script: ',name))
        source(file)
    }
    invisible()
}
