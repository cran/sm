".sm.Options" <-
    list(hmult = 1, h.weights = NA, period = NA,
         add = FALSE, band = NA, props = c(75, 50, 25), nbins = NA,
         positive = FALSE, delta = NA, display = NA,
         xlab = NA, ylab = NA, zlab = NA,
         xlim = NA, ylim = NA, zlim = NA, yht = NA,
         panel = FALSE, panel.plot = NA,
         ngrid = NA, eval.points = NA, rugplot = TRUE,
         col = NA, col.band = "cyan", col.mesh = "black",
         col.palette = topo.colors(12), col.points = "black",
         se = NA, se.breaks = c(-3, -2, 2, 3), lty = 1, pch = 1, cex = NA,
         theta = -30, phi = 40, size = 2, scaling = NULL,
         alpha = 0.7, alpha.mesh = 1, lit = FALSE,
         poly.index = 1, diff.ord = 2, test = NA, hull = TRUE, verbose = 1,
         df = NA, method = NA, structure.2d = "scaled", nboot = 100,
         describe = TRUE, show.script = TRUE, eval.grid = TRUE)

.onAttach <- function(library, pkg)
{
    ## we can't do this in .onLoad
    unlockBinding(".sm.Options", asNamespace("sm"))
    packageStartupMessage("Package `sm', version 2.2-4.1\n",
        "Copyright (C) 1997, 2000, 2005, 2007, 2008, A.W.Bowman & A.Azzalini\n",
        "Type help(sm) for summary information\n")
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
