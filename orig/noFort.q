.First.lib <- function(library, section)
{	
  assign(".sm.home", paste(library, section, sep="/"), w=0)
  cat("Library `sm'; Copyright (C) 1997 A.W.Bowman & A.Azzalini\n")
  cat("See `Licence' for details and conditions for use\n")
  invisible()
}

sm.density.3d <-
function(data, h=hnorm(data), hmult=1, contour = 75, ngrid = 20, 
	shadow = T, xlab = NA, ylab = NA, zlab = NA, theta = pi/4, phi = pi/4, 
	colour = T, title.colour = 3, label.colour = 1, axes.colour = 1, 
	plot.colour = 1, shadow.colour = 2, plt = T, cex = NA, maxpoly = 10000
	)
{
  stop("3D density estimation is not available in this version")
}
