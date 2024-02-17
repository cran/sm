sm <- function(x, y, data, subset, weights, bdeg = 3, h, model,
               random, ...) {
   
   data.missing    <- missing(data)
   subset.missing  <- missing(subset)
   weights.missing <- missing(weights)
   random.missing  <- missing(random)
   mf              <- if (data.missing) NULL else data
   
   if (!missing(y)) {
      x.name <- deparse(substitute(x))
      y.name <- deparse(substitute(y))
      if (weights.missing) weights <- NA
      if (missing(model)) model <- "none"
      return(sm.regression(x, y, h, model = model, weights = weights, ...))
   }
   else if (!inherits(x, "formula")) {
      x.name <- deparse(substitute(x))
      if (weights.missing) weights <- NA
      if (missing(model)) model <- "none"
      return(sm.density(x, h, model = model, weights = weights, xlab = x.name, ...))
   }
   
   # Insert redundant statements here to avoid NOTEs from devtools::check.
   reference <- NULL
   panel     <- NULL
   
   opt <- sm.options(list(...))
   replace.na(opt, display,   "none")
   replace.na(opt, reference, "none")
   replace.na(opt, panel,     FALSE)
   pam.formula    <- x
   
   if (!random.missing) {
      r.name <- deparse(substitute(random))
      random <- factor(as.vector(mf[ , r.name]))
   }
   
   #-----------------------------------------------------------------
   #      Parse model specification
   #-----------------------------------------------------------------
   
   terms.obj        <- terms(pam.formula, specials = "s")
   vars.inf         <- if (data.missing) eval.parent(attr(terms.obj, "variables"))
                       else eval(attr(terms.obj, "variables"), mf)
   term.labels      <- attr(terms.obj, "term.labels")
   s.ind            <- attr(terms.obj, "specials")$s
   response.ind     <- attr(terms.obj, "response")
   involved         <- attr(terms.obj, "factors")
   terms.linear     <- matrix(c(involved[s.ind, ]), ncol = length(term.labels))
   terms.linear     <- which(apply(terms.linear, 2, sum) == 0)
   nterms           <- length(term.labels)
   terms.smooth     <- which(!(1:nterms %in% terms.linear))
   rhs.linear       <- paste(term.labels[terms.linear], collapse = " + ")
   rhs.linear       <- if (nchar(rhs.linear) == 0) "1" else rhs.linear
   formula.linear   <- as.formula(paste("~", rhs.linear))
   names(vars.inf)  <- rownames(involved)
   bricks.type      <- sapply(vars.inf[-response.ind], mode)
   ind              <- (bricks.type == "numeric") &
      sapply(vars.inf[-response.ind], is.factor)
   bricks.type[ind] <- "factor"
   Xlinear          <- vars.inf[-response.ind][bricks.type != "list"]
   names(Xlinear)   <- names(bricks.type)[bricks.type != "list"]
   
   ylab             <- attr(terms.obj, "variables")
   ylab             <- strsplit(deparse(ylab), ",")[[1]][1]
   ylab             <- substr(ylab, 6, nchar(ylab))
   y                <- unlist(vars.inf[[response.ind]])
   X                <- list()
   xlabels          <- list()
   xlab             <- list()
   xdims            <- list()
   ndims            <- list()
   df               <- list()
   nseg             <- list()
   lambda           <- list()
   pord             <- list()
   period           <- list()
   increasing       <- list()
   xrange           <- list()
   fixed            <- list()
   fac              <- list()
   xmissing         <- FALSE
   
   if (length(terms.smooth) < 1) stop("there must be at least one smooth term.")
   
   if (any(apply(involved, 2, sum) > 3))
      stop("four-way interactions not yet implemented.")
   
   for (i in 1:length(terms.smooth)) {
      
      inv     <- which(involved[ , terms.smooth[i]] == 1)
      ilinear <- which(bricks.type[names(inv)] == "numeric")
      ifactor <- which(bricks.type[names(inv)] == "factor")
      if (length(ilinear) > 0)
         stop("interactions with linear terms are not yet implemented.")
      if (length(ifactor) > 1)
         stop("interactions with more than one factor are not yet implemented.")
      else if (length(ifactor) == 1) {
         fact     <- names(bricks.type)[ifactor]
         inv      <- inv[-match(fact, names(inv))]
         fac[[i]] <- Xlinear[[fact]]
      }
      else
         fac[[i]] <- NA
      
      nvars           <- length(inv)
      X[[i]]          <- matrix(nrow = length(y), ncol = 0)
      xlabels[[i]]    <- vector("character")
      xlab[[i]]       <- vector("character")
      xdims[[i]]      <- numeric()
      ndims[[i]]      <- numeric()
      df[[i]]         <- numeric()
      lambda[[i]]     <- numeric()
      pord[[i]]       <- numeric()
      period[[i]]     <- numeric()
      increasing[[i]] <- logical()
      nseg[[i]]       <- numeric()
      xrange[[i]]     <- matrix(nrow = 0, ncol = 2)
      fixed[[i]]      <- matrix(nrow = 0, ncol = 3)
      for (j in inv) {
         lambda[[i]]  <- c(lambda[[i]], vars.inf[[j]]$lambda)
         nseg[[i]]    <- c(nseg[[i]],   vars.inf[[j]]$nseg)
         xlabels[[i]] <- c(xlabels[[i]], vars.inf[[j]]$variables)
         newvar       <- if (data.missing) eval.parent(parse(text = vars.inf[[j]]$variables[1]))
         else              eval(parse(text = vars.inf[[j]]$variables[1]), mf)
         xdims[[i]]   <- if (is.matrix(newvar)) c(xdims[[i]], ncol(newvar)) else c(xdims[[i]], 1)
         if (length(vars.inf[[j]]$variables) > 1) {
            for (k in 2:length(vars.inf[[j]]$variables)) {
               newvar <- if (data.missing) cbind(newvar, eval.parent(parse(text = vars.inf[[j]]$variables[k])))
               else              cbind(newvar, eval(parse(text = vars.inf[[j]]$variables[k]), mf))
               xdims[[i]] <- c(xdims[[i]], ncol(newvar) - sum(xdims[[i]]))
            }
         }
         if (is.matrix(newvar)) {
            nms <- colnames(newvar)
            if (any(is.null(colnames(newvar))))
               nms <- paste(vars.inf[[j]]$variables,
                            "[", 1:ncol(newvar), "]", sep = "")
         }
         else
            nms <- vars.inf[[j]]$variables
         xlab[[i]]    <- c(xlab[[i]], nms)
         # xlab[[i]]    <- c(xlab[[i]], vars.inf[[j]]$variables)
         newvar       <- as.matrix(newvar)
         ndims.new    <- ncol(newvar)
         ndims[[i]]   <- c(ndims[[i]], ndims.new)
         pord[[i]]    <- c(pord[[i]], vars.inf[[j]]$pord)
         prd          <- vars.inf[[j]]$period
         if (length(prd) == 1 && is.na(prd)) prd <- rep(NA, ndims.new)
         if (length(prd) != ndims.new)
            stop("period does not match the columns of x.")
         period[[i]]  <- c(period[[i]], prd)
         if (any(!is.na(prd))) {
            for (k in 1:ndims.new)
               if (!is.na(prd[k])) newvar[ , k] <- newvar[ , k] %% prd[k]
         }
         incr            <- vars.inf[[j]]$increasing
         increasing[[i]] <- c(increasing[[i]], incr)
         xrng <- vars.inf[[j]]$xrange
         if ((ndims.new == 1) & (length(xrng) == 2))
            xrng <- matrix(xrng, nrow = 1)
         if (!is.matrix(xrng))
            xrng <- matrix(NA, nrow = ndims.new, ncol = 2)
         if (nrow(xrng) != ndims.new)
            stop("xrange does not match columns of x.")
         for (k in 1:ndims.new) {
            if (any(is.na(xrng[k, ]))) {
               if (!is.na(prd[k]))
                  xrng[k, ] <- c(0, prd[k])
               else
                  xrng[k, ] <- c(min(newvar[ , k], na.rm = TRUE), max(newvar[ , k], na.rm = TRUE))
               # xrange <- t(apply(xrange, 1, function(x) c(x[1] - 0.05 * diff(x), x[2] + 0.05 * diff(x))))
            }
         }
         xrange[[i]]  <- rbind(xrange[[i]], xrng)
         vinf         <- vars.inf[[j]]$fixed
         if (!is.matrix(vinf)) vinf <- matrix(vinf, nrow = 1)
         if (ncol(vinf) < 3) vinf <- cbind(vinf, 0)
         fixed[[i]]   <- rbind(fixed[[i]], vinf)
         X[[i]]       <- cbind(X[[i]], newvar)
         df.new       <- vars.inf[[j]]$df
         if (is.na(df.new)) df.new <- switch(ndims.new, 6, 12, 18)
         df[[i]]      <- c(df[[i]], df.new)
      }
      #      if (any(is.na(nseg[[i]])) | prod(nseg[[i]]) > 400)
      if (any(is.na(nseg[[i]])))
         nseg[[i]] <- rep(switch(sum(ndims[[i]]), 100, 17, 7), sum(ndims[[i]]))
      if (any(is.na(X[[i]]))) xmissing  <- TRUE
   }
   
   if (data.missing)
      B.linear <- model.matrix(formula.linear, parent.frame())
   else
      B.linear <- model.matrix(formula.linear, mf)
   if (nrow(B.linear) == 0)
      B.linear <- matrix(1, nrow = length(y), ncol = 1)

   if (opt$verbose > 1) {
      cat("Progress:\n")
      tim <- proc.time()
      timing <- function(tim) {
         elapsed <- round((proc.time() - tim)[1])
         unts <- "seconds"
         if (elapsed > 60) {
            elapsed <- round(elapsed / 60, 1)
            unts    <- if (elapsed > 1) "minutes" else "minute"
         }
         if (elapsed > 60) {
            elapsed <- round(elapsed / 60, 1)
            unts    <- "hours"
         }
         cat(elapsed, unts, "\n")
         proc.time()
      }
   }
   
   #-----------------------------------------------------------------
   #      Deal with subset and missing data
   #-----------------------------------------------------------------
   
   # Subset the data if required
   if (!subset.missing) {
      y <- y[subset]
      for (i in 1:length(X)) X[[i]] <- as.matrix(X[[i]][subset, ])
      B.linear <- as.matrix(B.linear[subset, ])
   }
   if (length(y) == 0)
      stop("no data after subsetting.")
   
   # Remove observations which have missing data
   ind.missing <- lapply(X, function(x) apply(x, 1, function(z) any(is.na(z))))
   ind.missing <- cbind(is.na(y), matrix(unlist(ind.missing), ncol = length(X)))
   ind.missing <- cbind(ind.missing, apply(B.linear, 1, function(z) any(is.na(z))))
   if (!random.missing)
      ind.missing <- cbind(ind.missing, is.na(random))
   ind.missing <- apply(ind.missing, 1, any)
   if (all(ind.missing))
      stop("no data after the removal of missing values.")
   if (any(ind.missing)) {
      y <- y[!ind.missing]
      for (i in 1:length(X)) X[[i]] <- as.matrix(X[[i]][!ind.missing, ])
      B.linear <- B.linear[!ind.missing, , drop = FALSE]
      if (!random.missing)
         random   <- random[!ind.missing]
      if (opt$verbose > 0) cat("warning: missing data removed.\n")
   }
   
   #-----------------------------------------------------------------
   #      Construct matrices
   #-----------------------------------------------------------------
   
   if (opt$verbose > 1) cat(" constructing matrices .......... ")
   
   P <- list(length = length(terms.smooth))
   B <- B.linear
   m <- ncol(B)

   for (i in 1:length(terms.smooth)) {
      mat    <- ps.matrices(X[[i]], xrange[[i]], ndims = ndims[[i]],
                            nseg = nseg[[i]], pord = pord[[i]], period = period[[i]])
      if (all(is.na(fac[[i]]))) {
         B      <- cbind(B, mat$B)
         m      <- c(m, ncol(mat$B))
         P[[i]] <- mat$P
      }
      else {
         Btemp <- matrix(nrow = length(y), ncol = 0)
         for (j in levels(fac[[i]]))
            Btemp <- cbind(Btemp, mat$B * as.numeric(fac[[i]] == j))
         B      <- cbind(B, Btemp)
         m      <- c(m, ncol(Btemp))
         nlevs  <- length(levels(fac[[i]]))
         pdim   <- nlevs * ncol(mat$P)
         P[[i]] <- matrix(0, nrow = pdim, ncol = pdim)
         for (j in 1:nlevs) {
            ind <- (j - 1) * ncol(mat$P) + (1:ncol(mat$P))
            P[[i]][ind, ind] <- mat$P
         }
         P[[i]] <- P[[i]] + matrix(1, ncol = nlevs, nrow = nlevs) %x% diag(ncol(mat$B))
      }
      xrange[[i]] <- mat$xrange
   }
   
   if (opt$verbose > 1) tim <- timing(tim)
   
   #-----------------------------------------------------------------
   #      Construct matrix products
   #-----------------------------------------------------------------
   
   if (opt$verbose > 1) cat(" constructing matrix products ... ")
   
   b.ind <- list(length = length(m))
   for (i in 1:length(terms.smooth))
      b.ind[[i]] <- (cumsum(m)[i] + 1):cumsum(m)[i + 1]
   
   if (weights.missing) {
      btb   <- crossprod(B)
      bty   <- t(B)  %*% y
   }
   else if (is.vector(weights)) {
      btb   <- t(B * weights)  %*% B
      bty   <- t(B * weights)  %*% y
   }
   else if (is.matrix(weights)) {
      btb   <- t(B) %*% weights %*% B
      bty   <- t(B) %*% weights %*% y
   }
   else
      stop("the weights argument is inappropriate.")
   
   if (opt$verbose > 1) tim <- timing(tim)
   
   #-----------------------------------------------------------------
   #      Select the smoothing parameters (if required)
   #-----------------------------------------------------------------
   
   if (opt$verbose > 1) cat(" setting smoothness ............. ")
   
   for (i in 1:length(terms.smooth)) {
      # code doesn't currently handle more than one df for terms with more than one variable.
      df[[i]] <- sum(df[[i]])
      if (df[[i]] > prod(nseg[[i]] + 3))
         stop(paste("df is too large for the value of nseg in term", i))
      if (any(is.na(lambda[[i]])))
         lambda[[i]] <- lambda.select(btb[b.ind[[i]], b.ind[[i]]], P[[i]], df[[i]])
   }
   
   if (opt$verbose > 1) tim <- timing(tim)
   
   #-----------------------------------------------------------------
   #      Fit the model
   #-----------------------------------------------------------------
   
   if (opt$verbose > 1) cat(" fitting ........................ ")
   
   Pall  <- matrix(0, nrow = ncol(B), ncol = ncol(B))
   for (i in 1:length(terms.smooth))
      Pall[b.ind[[i]], b.ind[[i]]] <- lambda[[i]] * P[[i]]
   B1       <- solve(btb + Pall)
   alpha    <- as.vector(B1 %*% bty)
   df.model <- sum(btb * t(B1))
   df.error <- length(y) - sum(btb * (2 * B1 - B1 %*% btb %*% B1))
   
   if (opt$verbose > 1) tim <- timing(tim)
   
   #-----------------------------------------------------------------
   #      Constraints (if requested)
   #-----------------------------------------------------------------
   
   # Force the estimate to pass through fixed points
   # if (length(terms.smooth) == 1 & ndims[[1]] == 1 & all(!is.na(fixed[[1]]))) {
   if (length(terms.smooth) == 1 & ndims[[1]] <= 2 & all(!is.na(fixed[[1]]))) {
      fxd <- fixed[[1]]
      ind <- any(fxd[ , 1] < xrange[[1]][1, 1] - 100*.Machine$double.eps) |
         any(fxd[ , 1] > xrange[[1]][1, 2] + 100*.Machine$double.eps)
      if (ndims[[1]] == 2)
         ind <- ind | any(fxd[ , 2] < xrange[[1]][2, 1] - 100*.Machine$double.eps) |
         any(fxd[ , 2] > xrange[[1]][2, 2] + 100*.Machine$double.eps)
      if (ind) stop("fixed points must be inside the range of the data.")
      fx <- fxd[ , 1:ndims[[1]]]
      fx <- matrix(c(fx), ncol = ndims[[1]])
      A     <- cbind(1, ps.matrices(fx, xrange[[1]], ndims[[1]], nseg[[1]])$B,
                     pord = pord[[1]])
      alpha <- alpha +  B1 %*% t(A) %*% solve(A %*% B1 %*% t(A)) %*%
         (fxd[ , ndims[[1]] + 1] - A %*% alpha)
   }
   
   # Increasing function: specific to nseg 17 and so 20 basis fns in each dim.
   if (random.missing && (length(ndims) == 1) && ndims[[1]] == 2 &&
       increasing[[1]] && all(nseg[[1]] == 17) && bdeg == 3) {
      D1    <- diff(diag(20)) %x% diag(20)
      D2    <- diag(20)       %x% diff(diag(20))
      delta <- 1
      while (delta > 1e-5) {
         v1        <- as.numeric(c(D1 %*% alpha[-1]) <= 0)
         v2        <- as.numeric(c(D2 %*% alpha[-1]) <= 0)
         mat1      <- 100 * lambda[[1]] * t(D1) %*% diag(v1) %*% D1
         mat2      <- 100 * lambda[[1]] * t(D2) %*% diag(v2) %*% D2
         mat       <- matrix(0, nrow = ncol(B), ncol = ncol(B))
         mat[2:401, 2:401] <- mat1 + mat2
         B1        <- solve(btb + Pall + mat)
         alpha.old <- alpha
         alpha     <- as.vector(B1 %*% bty)
         delta     <- sum((alpha - alpha.old)^2) / sum(alpha.old^2)
      }
      mu       <- c(B %*% alpha)
      df.model <- NULL
      df.error <- NULL
   }
   
   #-----------------------------------------------------------------
   #      Summaries (and random effects if requested)
   #-----------------------------------------------------------------
   
   if (!random.missing) {
      if (opt$verbose > 1) cat(" estimating the random effect ... ")
      nrnd   <- nlevels(random)
      n      <- length(y)
      ind    <- cbind(1:n, as.numeric(random))
      utu    <- diag(table(as.numeric(random)))
      btu    <- t(apply(B, 2, function(x) tapply(x, random, sum)))
      uty    <- tapply(y, random, sum)
      sig    <- 1
      sigr   <- 1
      sigo   <- 2
      sigro  <- 2
      p1     <- alpha
      p2     <- rep(0, nrnd)
      m1     <- c(B1 %*% bty)
      M1     <- B1 %*% btu
      while (abs(sig - sigo)/sigo > 1e-6 | abs(sigr - sigro)/sigro > 1e-6) {
         invu  <- solve(utu + diag(rep(sig^2/sigr^2, nrnd)))
         m2    <- as.vector(invu %*% uty)
         M2    <- invu %*% t(btu)
         p1    <- solve(diag(sum(m)) - M1 %*% M2) %*% (m1 - M1 %*% m2)
         p2    <- solve(diag(nrnd)   - M2 %*% M1) %*% (m2 - M2 %*% m1)
         ress  <- sum(y^2) + c(t(p1) %*% btb %*% p1 + t(p2) %*% utu %*% p2 - 2 * t(p1) %*% bty +
                                  2 * t(p1) %*% btu %*% p2 - 2 * t(p2) %*% uty)
         nu    <- sum(diag(invu)) / sigr^2
         sigo  <- sig
         sigro <- sigr
         sig   <- sqrt(ress / (n - nrnd + nu))
         sigr  <- sqrt(sum(p2^2) / (nrnd - nu))
         # print(c(sig, sigr))
      }
      alpha     <- p1
      sigma     <- sig
      sigma.r   <- sigr
      N1        <- solve(diag(sum(m)) - M1 %*% M2)
      U1        <- solve(utu + diag(nrnd) * sigma^2 / sigma.r^2)
      cov.alpha <- sigma^2 * btb -
         sigma^2 * btu %*% U1 %*% t(btu) +
         sigma.r^2 * btu %*% t(btu) -
         sigma.r^2 * btu %*% utu %*% U1 %*% t(btu) -
         sigma^2 * btu %*% U1 %*% t(btu) +
         sigma^2 * btu %*% U1 %*% utu %*% U1 %*% t(btu) -
         sigma.r^2 * btu %*% U1 %*% utu %*% t(btu) +
         sigma.r^2 * btu %*% U1 %*% utu %*% utu %*% U1 %*% t(btu)
      cov.alpha <- N1 %*% B1 %*% cov.alpha %*% B1 %*% t(N1)
      rss       <- NULL
      R.squared <- NULL
      mu         <- c(B %*% alpha)
      if (opt$verbose > 1) tim <- timing(tim)
   }
   else {
      mu         <- c(B %*% alpha)
      sigma      <- sqrt(sum((y - mu)^2) / df.error)
      cov.alpha  <- B1 %*% btb %*% t(B1) * sigma^2
      rss        <- sum((y - mu)^2)
      tss        <- sum((y - mean(y))^2)
      R.squared  <- 100 * (tss - rss) / tss
   }
   
   # P[[1]] <- P[[1]] %x% diag(m[2])
   # P[[2]] <- diag(m[2]) %x% P[[2]]
   
   # If there is only one term, include the mean
   # if (nterms == 1) b.ind[[1]] <- c(1, b.ind[[1]])
   
   result <- list(fitted = mu, alpha = alpha, m = m, B = B,
                  bty = bty, btb = btb, B1 = B1, Pall = Pall, xlabels = xlabels,
                  linear.matrix = B.linear, n = nrow(B),
                  terms.linear = terms.linear, terms.smooth = terms.smooth,
                  xlab = xlab, ylab = ylab, term.labels = term.labels,
                  lambda = lambda, ndims = ndims, xdims = xdims,
                  y = y, X = X, fac = fac, Xlinear = Xlinear,
                  bricks.type = bricks.type,
                  sigma = sigma, cov.alpha = cov.alpha, b.ind = b.ind,
                  df = df, df.model = df.model, df.error = df.error,
                  rss = rss, R.squared = R.squared, xrange = xrange,
                  nseg = nseg, bdeg = bdeg, pord = pord, period = period,
                  increasing = increasing, pam.formula = pam.formula,
                  formula.linear = formula.linear,
                  involved = involved, nterms = nterms,
                  ind.missing = which(ind.missing))
   if (!weights.missing) result$weights <- weights
   if (!random.missing)  result$sigma.r <- sigma.r
   class(result) <- "pam"
   
   # if (nterms == 1 & ndims[[1]] <= 2) {
   # 	  if (missing(model)) model <- "none"
   # 	  colnames(result$X[[1]]) <- xlab[[1]]
   #    if (ndims[[1]] == 1) sm.regression(result$X[[1]], result$y, h = h, model = model,
   #                                          xlab = xlab[[1]], ylab = ylab, ...)
   #    else                 sm.regression(result$X[[1]], result$y, h = h, model = model, zlab = ylab, ...)
   #    return(invisible(result))
   # }
   # else
   if (opt$display != "none") plot(result, ...)
   
   # if (nterms == 1 & ndims[[1]] <= 2) {
   # if (opt$panel) {
   # replace.na(opt, df.max, switch(ndims[[1]], 20, 50, 100))
   # df.min  <- switch(ndims[[1]], 2, 4, 8) + 0.1
   # df.max  <- if (!opt$panel) df[[1]] else min(length(y) - 5, opt$df.max)
   # df.min  <- if (!opt$panel) df[[1]] else df.min
   # Pall    <- rbind(0, cbind(0, P[[1]]))
   # llambda    <- 0
   # llambda.df <- lambda.df(exp(max(llambda)), btb, Pall)
   # while (min(llambda.df) >= df.min) {
   # llambda    <- c(llambda, max(llambda) + log(10))
   # llambda.df <- c(llambda.df, lambda.df(exp(max(llambda)), btb, Pall))
   # }
   # while (max(llambda.df) <= df.max) {
   # llambda    <- c(llambda, min(llambda) - log(10))
   # llambda.df <- c(llambda.df, lambda.df(exp(min(llambda)), btb, Pall))
   # }
   # df.fun <- approxfun(llambda.df, llambda)
   
   # sm.pam.draw <- function(pam.panel) {
   # plot(pam.panel$model, options = pam.panel$opt)
   # # title(pam.panel$df)
   # pam.panel
   # }
   # sm.pam.redraw <- function(pam.panel) {
   # # pam.panel$model$lambda <- lambda.select(pam.panel$model$btb, pam.panel$model$bty,
   # #                                         Pall, pam.panel$df)
   # pam.panel$model$lambda <- exp(pam.panel$df.fun(pam.panel$df))
   # B1 <- solve(pam.panel$model$btb + pam.panel$model$lambda * pam.panel$Pall)
   # pam.panel$model$alpha  <- as.vector(B1 %*% pam.panel$model$bty)
   # pam.panel$model$fitted <- c(pam.panel$model$B %*% pam.panel$model$alpha)
   # pam.panel$opt$se       <- pam.panel$se
   # pam.panel$opt$theta    <- pam.panel$theta
   # pam.panel$opt$phi      <- pam.panel$phi
   # rp.control.put(pam.panel$panelname, pam.panel)
   # rp.tkrreplot(pam.panel, plot)
   # pam.panel
   # }
   # opt1 <- opt
   # opt1$panel <- FALSE
   # pam.panel <- rp.control(model = result, opt = opt1, Pall = rbind(0, cbind(0, P[[1]])),
   # df = opt$df, df.fun = df.fun, theta = opt$theta, phi = opt$phi)
   # rp.tkrplot(pam.panel, plot, sm.pam.draw, hscale = opt$hscale, vscale = opt$vscale, pos = "right")
   # rp.slider(pam.panel, df, df.min, df.max, sm.pam.redraw, showvalue = TRUE)
   # rp.checkbox(pam.panel, se, sm.pam.redraw, title = "Standard errors")
   # if (ndims[[1]] == 2) {
   # rp.slider(pam.panel, theta, -180, 180, sm.pam.redraw, "persp angle 1")
   # rp.slider(pam.panel, phi,      0,  90, sm.pam.redraw, "persp angle 2")
   # }
   # }
   # else if (opt$display != "none")
   # plot(result, ...)
   # }
   
   invisible(result)
}

#----------------------------------------------------------------------------
#                             lambda.select
#----------------------------------------------------------------------------

lambda.df <- function(lambda, btb, P, df) {
   B1   <- solve(btb + lambda * P)
   sum(diag(btb %*% B1)) - df
}

lambda.select <- function(btb, P, df, method = "df") {
   #    This currently uses the same lambda in all dimensions
   if (method == "df") {
      lambda <- c(1, 1)
      f      <- rep(lambda.df(lambda[1], btb, P, df), 2)
      mult   <- if (f[1] < 0) 0.1 else 10
      while (sign(f[1]) == sign(f[2])) {
         lambda <- c(lambda[2], lambda[2] * mult)
         f      <- c(f[2], lambda.df(lambda[2], btb, P, df))
      }
      if (diff(f) > 0) {
         lambda <- rev(lambda)
         f      <- rev(f)
      }
      lambda <- uniroot(lambda.df, interval = lambda, btb, P, df,
                        f.lower = f[1], f.upper = f[2])$root
      
      # 	  lambda <- 1
      # 	  f.new  <- 0
      # 	  while (f.new <= df) {
      # 	  	lower   <- lambda
      # 	  	f.lower <- f.new
      # 	  	lambda  <- lambda / 10
      # 	  	f.new   <- lambda.df(lambda, btb, P)
      # 	  }
      # 	  lambda <- 1
      # 	  f.new  <- df + 1
      # 	  while (f.new >= df) {
      # 	  	upper   <- lambda
      # 	  	f.upper <- f.new
      # 	  	lambda  <- lambda * 10
      # 	  	f.new   <- lambda.df(lambda, btb, P)
      # 	  }
      # print(c(lower, upper))
      # #
      # # 	  lower <- lambda / 10
      # #  	lambda  <- 0.1
      # #  	f.upper <- df + 1
      # #  	while (f.upper >= df) {
      # #  		lambda  <- lambda * 10
      # #  		f.upper <- lambda.df(lambda, btb, P)
      # #  	}
      # # 	 upper <- lambda * 10
      # # 	 f.lower <- lambda.df(lambda, btb, P)
      # #    lower  <- lambda
      # #    lambda <- 1
      # #    f.lambda <- df - 1
      # #    while (lambda.df(lambda, btb, P) >= df) lambda <- lambda * 10
      # #    upper  <- lambda
      #    lambda.crit <- function(lambda, btb, P, df)
      #       lambda.df(lambda, btb, P) - df
      #    cat("entering uniroot ... \n")
      #    result <- uniroot(lambda.crit, interval = c(lower, upper),
      #    									f.lower = f.lower, f.upper = f.upper, btb, P, df)
      #    # cat("result$root", result$root, "\n")
      #    lambda <- result$root
      #    cat("leaving uniroot ... \n")
   }
   lambda
}

#----------------------------------------------------------------------------
#                                  s
#----------------------------------------------------------------------------

s <- function(..., lambda = NA, df = NA, period = NA, increasing = FALSE,
              xrange = NA, nseg = NA, pord = 2, fixed = c(NA, NA, NA)) {
   vars.list <- as.list(substitute(list(...)))[-1]
   nvar <- length(vars.list)
   if (nvar > 3)
      stop("smooth terms can be constructed from only 1, 2 or 3 variables.")
   variables <- character(0)
   for (i in 1:nvar) variables <- c(variables, deparse(vars.list[[i]]))
   list(variables = variables, lambda = lambda, df = df, period = period,
        increasing = increasing, xrange = xrange, nseg = nseg, pord = pord, fixed = fixed)
}

#----------------------------------------------------------------------------
#                                  ps.matrices
#----------------------------------------------------------------------------

ps.matrices <- function(x, xrange, ndims, nseg, bdeg = 3, pord = 2, period = NA,
                        decompose =  TRUE, penalty = TRUE) {
   
   # Compute a set of basis functions and a penalty matrix associated with x.
   # An intercept term and the main effects of any interaction terms are removed.
   
   n     <- nrow(x)
   ndimx <- ncol(x)
   if (ndimx > 3) stop("terms with more than three dimensions cannot be used.")
   
   if (missing(nseg)) nseg <- rep(switch(ndimx, 100, 17, 7), ndimx)
   if (length(bdeg) < ndimx) bdeg <- rep(bdeg[1], ndimx)
   
   # Compute B-spline basis
   
   b <- list(length = ndimx)
   m <- vector(length = ndimx)
   for (i in 1:ndimx) {
      b[[i]] <- bbase(x[,i], xl = xrange[i , 1], xr = xrange[i, 2], nseg = nseg[i],
                      deg = bdeg[i])
      m[i]   <- ncol(b[[i]])
   }
   
   B <- b[[1]]
   if (ndimx > 1)
      B <- t(apply(cbind(b[[1]], b[[2]]), 1,
                   function(x) c(x[1:m[1]] %x% x[-(1:m[1])])))
   if (ndimx == 3)
      B <- t(apply(cbind(B,  b[[3]]), 1,
                   function(x) c(x[1:(m[1]*m[2])] %x% x[-(1:(m[1]*m[2]))])))
   
   result <- list(B = B, xrange = xrange, nseg = nseg, bdeg = bdeg, pord = pord)
   if (!penalty) return(invisible(result))
   
   # Construct smoothness penalty matrices
   P <- list()
   for (i in 1:ndimx) {
      for (j in 1:length(pord)) {
         Pi <- diff(diag(m[i]), diff = pord[j])
         if (!is.na(period[i])) {
            z  <- c(1, rep(0, m[i] - 4), -1)
            Pi <- rbind(Pi, c(z, 0, 0))
            Pi <- rbind(Pi, c(0, z, 0))
            Pi <- rbind(Pi, c(0, 0, z))
         }
         Pi <- crossprod(Pi)
         if (j == 1) P[[i]] <- Pi
         else        P[[i]] <- P[[i]] + Pi
      }
   }
   if (ndimx >= 2) {
      P[[1]] <- P[[1]] %x% diag(m[2])
      P[[2]] <- diag(m[1]) %x% P[[2]]
   }
   if (ndimx == 3) {
      P[[1]] <- P[[1]] %x% diag(m[3])
      P[[2]] <- P[[2]] %x% diag(m[3])
      P[[3]] <- diag(m[1]) %x% diag(m[2]) %x% P[[3]]
   }
   pmat <- matrix(0, nrow = ncol(B), ncol = ncol(B))
   for (i in 1:ndimx)
      pmat <- pmat + P[[i]]
   
   #     Construct anova constraint penalty matrices
   if (length(ndims) == 1) {
      # Sum of coefficients constraint
      # cmat <- matrix(1, nrow = prod(m), ncol = prod(m))
      # Sum of estimated values constraint
      Bsum <- apply(B, 2, sum)
      cmat <- Bsum %o% Bsum
      # Corner point constraint (first coefficient is 0
      # cmat <- diag(c(1, rep(0, ncol(B) - 1)))
      pmat <- pmat + cmat
   }
   else if (length(ndims) == 2) {
      if (all(ndims == c(1, 1))) ind <- c(m[1], m[2])
      if (all(ndims == c(1, 2))) ind <- c(m[1], m[2] * m[3])
      if (all(ndims == c(2, 1))) ind <- c(m[1] * m[2], m[3])
      pmat <- pmat + matrix(1, nrow = ind[1], ncol = ind[1]) %x% diag(ind[2])
      pmat <- pmat + diag(ind[1]) %x% matrix(1, nrow = ind[2], ncol = ind[2])
   }
   else if (length(ndims) == 3) {
      pmat <- pmat + matrix(1, nrow = m[1], ncol = m[1]) %x% diag(m[2]) %x% diag(m[3])
      pmat <- pmat + diag(m[1]) %x% matrix(1, nrow = m[2], ncol = m[2]) %x% diag(m[3])
      pmat <- pmat + diag(m[1]) %x% diag(m[2]) %x% matrix(1, nrow = m[3], ncol = m[3])
   }
   
   result$P <- pmat
   if (length(ndims) == 1) result$cmat <- cmat
   invisible(result)
}

bbase <- function(x, xl = min(x), xr = max(x), nseg = 10, deg = 3) {
   # Construct B-spline basis
   dx <- (xr - xl) / nseg
   knots <- seq(xl - deg * dx, xr + deg * dx, by = dx)
   P <- outer(x, knots, tpower, deg)
   n <- dim(P)[2]
   D <- diff(diag(n), diff = deg + 1) / (gamma(deg + 1) * dx ^ deg)
   B <- (-1) ^ (deg + 1) * P %*% t(D)
   B
}

tpower <- function(x, t, p)
   # Truncated p-th power function
   (x - t) ^ p * (x > t)
