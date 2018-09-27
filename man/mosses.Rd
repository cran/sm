\name{mosses}

\alias{mosses}

\title{Heavy metals in mosses in Galicia.}

\description{Mosses are used as a means of measuring levels of heavy metal concentrations in the atmosphere, since most of the nutrient uptake of the mosses is from the air.  This technique for large-scale monitoring of long-range transport processes has been used in Galicia, in North-West Spain, over the last decade, as described by Fernandez et al. (2005).  In 2006, in both March and September, measurements of different metals were collected at 148 points lying almost in a regular grid over the region with 15 km spacing in north-south and east-west directions.  According to the ecologists' expertise, the period between the two samples, passing from a humid to a dry season, is enough time to guarantee the independence of the observed processes.

The dataset consists of a list with six components
\tabular{ll}{
\code{loc.m} \tab a two-column matrix containing grid locations of the March monitoring sites \cr
\code{loc.s} \tab a two-column matrix containing grid locations of the September monitoring sites \cr
\code{Co.m}  \tab cobalt concentration (log scale) in March \cr
\code{Co.s}  \tab cobalt concentration (log scale) in September \cr
\code{Hg.m}  \tab mercury concentration (log scale) in March \cr
\code{Hg.s}  \tab mercury concentration (log scale) in September \cr}

Source: The data were kindly made available by the Ecotoxicology and Vegetal Ecophysiology research group in the University of Santiago de Compostela.
}

\references{Fernandez, J., Real, C., Couto, J., Aboal, J., Carballeira, A. (2005).
The effect of sampling design on extensive biomonitoring surveys of air pollution.
Science of the Total Environment, 337, 11-21.}

\examples{
\dontrun{
# Comparison of Co in March and September
   
with(mosses, {
	
   nbins <- 12
   vgm.m <- sm.variogram(loc.m, Co.m, nbins = nbins, original.scale = TRUE,
                        ylim = c(0, 1.5))
   vgm.s <- sm.variogram(loc.s, Co.s, nbins = nbins, original.scale = TRUE,
                        add = TRUE, col.points = "blue")
                        
   trns <- function(x) (x / 0.977741)^4
   del <- 1000
   plot(vgm.m$distance.mean, trns(vgm.m$sqrtdiff.mean), type = "b",
         ylim = c(0, 1.5), xlab = "Distance", ylab = "Semi-variogram")
   points(vgm.s$distance.mean - del, trns(vgm.s$sqrtdiff.mean), type = "b",
         col = "blue", pch = 2, lty = 2)

   plot(vgm.m$distance.mean, trns(vgm.m$sqrtdiff.mean), type = "b",
         ylim = c(0, 1.5), xlab = "Distance", ylab = "Semi-variogram")
   points(vgm.s$distance.mean - del, trns(vgm.s$sqrtdiff.mean), type = "b",
         col = "blue", pch = 2, lty = 2)
   segments(vgm.m$distance.mean, trns(vgm.m$sqrtdiff.mean - 2 * vgm.m$se),
         vgm.m$distance.mean, trns(vgm.m$sqrtdiff.mean + 2 * vgm.m$se))
   segments(vgm.s$distance.mean - del, trns(vgm.s$sqrtdiff.mean - 2 * vgm.s$se),
         vgm.s$distance.mean - del, trns(vgm.s$sqrtdiff.mean + 2 * vgm.s$se),
         col = "blue", lty = 2)

   mn <- (vgm.m$sqrtdiff.mean + vgm.s$sqrtdiff.mean) / 2
   se <- sqrt(vgm.m$se^2 + vgm.s$se^2)
   plot(vgm.m$distance.mean, trns(vgm.m$sqrtdiff.mean), type = "n",
        ylim = c(0, 1.5), xlab = "Distance", ylab = "Semi-variogram")
   polygon(c(vgm.m$distance.mean, rev(vgm.m$distance.mean)),
        c(trns(mn - se), rev(trns(mn + se))),
        border = NA, col = "lightblue")  
   points(vgm.m$distance.mean, trns(vgm.m$sqrtdiff.mean))
   points(vgm.s$distance.mean, trns(vgm.s$sqrtdiff.mean), col = "blue", pch = 2)

   vgm1 <- sm.variogram(loc.m, Co.m, nbins = nbins, varmat = TRUE, 
                        display = "none")
   vgm2 <- sm.variogram(loc.s, Co.s, nbins = nbins, varmat = TRUE,
                        display = "none")

   nbin  <- length(vgm1$distance.mean)
   vdiff <- vgm1$sqrtdiff.mean - vgm2$sqrtdiff.mean
   tstat <- c(vdiff \%*\% solve(vgm1$V + vgm2$V) \%*\% vdiff)
   pval  <- 1 - pchisq(tstat, nbin)
   print(pval)
})

# Assessing isotropy for Hg in March

with(mosses, {
   sm.variogram(loc.m, Hg.m, model = "isotropic")
})

# Assessing stationarity for Hg in September

with(mosses, {
   vgm.sty <- sm.variogram(loc.s, Hg.s, model = "stationary")
   i <- 1
   image(vgm.sty$eval.points[[1]], vgm.sty$eval.points[[2]], vgm.sty$estimate[ , , i],
         col = topo.colors(20))
   contour(vgm.sty$eval.points[[1]], vgm.sty$eval.points[[2]], vgm.sty$sdiff[ , , i],
         col = "red", add = TRUE)
})

}
}

\keyword{smooth}
\keyword{regression}
