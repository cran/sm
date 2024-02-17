".sm.Options" <-
    list(hmult = 1, h.weights = NA, period = NA,
         add = FALSE, band = NA, props = c(75, 50, 25), nbins = NA,
         positive = FALSE, delta = NA, display = NA, 
         hscale = 1, vscale = 1,
         xlab = NA, ylab = NA, zlab = NA, x1lim = NA, x2lim = NA,
         xlim = NA, ylim = NA, zlim = NA, yht = NA, asp = NA,
         cex = NA, cex.axis = NA, cex.lab = NA, labcex = NA, key = TRUE,
         model = "none", reference = "none",
         panel = FALSE, panel.plot = NA,
         ngrid = NA, eval.points = NA, rugplot = TRUE, 
         col = NA, col.band = "cyan", col.mesh = "black", 
         col.points = "black",
         col.palette = topo.colors(12), col.palette.fn = topo.colors, 
         superimpose = NA,
         se = NA, se.breaks = c(-3, -2, 2, 3), lty = 1, lwd = 1, pch = 1,
         theta = -30, phi = 40, size = 2, scaling = NULL, 
         alpha = 0.7, alpha.mesh = 1, lit = FALSE,
         poly.index = 1, diff.ord = 2, test = NA, hull = TRUE, verbose = 1, 
         df = NA, df.max = NA, method = NA, structure.2d = "scaled", nboot = 100, 
         describe = TRUE, show.script = TRUE, eval.grid = TRUE,
         mask.method = "hull", partial.residuals = TRUE, nlevels = 20,
         include.mean = NA, include.terms = "single", include.lower.reference = FALSE,
         transform.response = I, eqscplot = FALSE, order = 1:3,
         colour.key = NA, z.key = NA,
         deriv = NA, deriv.order = NA)


.onAttach <- function(library, pkg)
{
    ## we can't do this in .onLoad
    unlockBinding(".sm.Options", asNamespace("sm"))
    packageStartupMessage("Package 'sm', version 2.2-6.0: type help(sm) for summary information")
    invisible()
}


isMatrix <- function(x) length(dim(x)) == 2

isInteger <- function(x) all(x == round(x))

sm.script <- function(name, path) {
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
