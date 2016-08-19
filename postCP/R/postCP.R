postcp <-
function(formula, data, bp=integer(), family = gaussian(), sigma=1.0, maxFB=FALSE)
{
  #Extract data from environment if data is missing
  if (missing(data))
    data <- environment(formula)
  X <- model.matrix(formula, data)          #covariate matrix
  n <- nrow(X)                              #length of the data
  p <- ncol(X)                              #number of variables
  cp <- bp                                  #change points
  CP <- c(0, cp, n)                         #change points vector formatted
  K <- length(bp) + 1                       #number of segments
  beta <- matrix(0, p, K)                   #matrix of initial regression coefficients
  y <- rep(NA, n)
  mf <- model.frame(formula, data)
  y <- model.response(mf, "any")            #response variable

  #Response variable not specified
  if (any(is.null(y))) {
    stop("Response variable not specified in the formula");
  }

  #Error handling checks for change points begin here ----------------------------
  table.seg <- table(bp)
  if (sum(table.seg > 1) > 0) {
    stop(paste("Change-point",names(table.seg)[table.seg > 1],"is repeated at least twice, please include only once\n"));
  }
  if (n < 2){
    stop("Not enough data to segment");
  }
  if (K > n) {
    stop("Number of segments larger than data");
  }
  if (length(bp) > 0) {
    if (max(bp) >= n) stop("At least one change-point initialized to larger than n");
    if (min(bp) < 1)    stop("Location of first change-point must be at least 1");
    if (sum(round(bp) != bp) > 0) stop("Enter all change-points as integers");
  }
  #Error handling checks for change points end here -------------------------------

  #calculation of the initial estimate for beta
  for (k in 1:K) {
      coeffK <- glm(formula=formula, data=data[(CP[k] + 1):(CP[k+1]), ], family=family);
      beta[, k] <- coefficients(coeffK)
  }

  #Check for linear dependency of variables
  if (any(is.na(beta))) {
    na_rows <- which(is.na(beta), arr.ind=TRUE)[, 1]
    stop(paste("\"", colnames(X)[na_rows], "\"", "variable specified in the formula is linearly dependent\n"))
  }

  if (nrow(beta)!= p) {
    stop("Dimensions of the data matrix do not match with parameters(beta matrix)")
  }
  Xbeta <- X %*% beta

  # compute response variable based on the family
  if (is.character(family)) {
    family <- get(family, mode = "function", envir = parent.frame())
  }
  if (is.function(family)) {
    family <- family()
  }

  if (family$family == 'poisson') {
    if (sum(Xbeta < 0) > 0) {
      stop("Negative numbers with Poisson distribution specified, choose family=gaussian");
    }
    if (sum(round(Xbeta) != Xbeta)>0) {
      stop("Non-integer numbers with Poisson distribution specified, choose family=gaussian");
    }

    # compute the log-evidence le
    # this is basically the only model-specific part !
    le <- y * log(Xbeta) - y - lgamma(y+1)
  } else if (family$family == 'gaussian') {
    # compute the log-evidence le
    # this is basically the only model-specific part !
    le <- -(Xbeta - y)^2/2/sigma^2 - log(sigma) - 0.5 * log(2 * pi)
  } else if (family$family == 'binomial') {
    if(sum(Xbeta < 0) > 0) {
      stop("Negative numbers with Binomial distribution specified, choose family=gaussian");
    }

    # compute the log-evidence le
    # this is basically the only model-specific part !
    le <-  y * log(Xbeta) + (1 - y) * log(1 - Xbeta)
  } else if (family$family == 'Gamma') {
    if (sum(Xbeta < 0) > 0) {
      stop("Negative numbers with Gamma distribution specified, choose family=gaussian");
    }
    # compute the log-evidence le
    # this is basically the only model-specific part !
    le <- (Xbeta[1] - 1) * log(y) - y/Xbeta[2] - Xbeta[1] * log(Xbeta[2]) - lgamma(Xbeta[1])

  } else {
    print(family)
    stop("'family' not recognized")
  }

  # build workspace for forward-backward
  lFw <- matrix(0, n, K)
  lFw[] <- le[]
  lBk <- matrix(0, n, K)
  # call forward-backward
  if (maxFB) {  # call maxFwBk if the most probable change point configuration is required
    loglik <- maxFwBk(lFw, lBk) - lchoose(n-1, K-1)
  } else { # call FwBk if marginal distributions for change points are required
    loglik <- FwBk(lFw, lBk) - lchoose(n-1, K-1)
  }
  # compute posterior distribution
  post <- lFw + lBk
  post <- exp(post - apply(post, 1, logsumexp))
  # compute posterior cp distribution
  post.cp <- matrix(0, n, K-1)
  for (k in 1:(K-1)) {
    aux <- lFw[-n, k] + le[-1, k+1] + lBk[-1, k+1]
    post.cp[-n, k] <- exp(aux-max(aux))
  }

  # parameter estimation
  le.updated <- matrix(0, n, K)
  param.updated <- matrix(0, p, K)
  for (k in 1:(K)) {
    fit <- glm(y ~ X-1, weights=post[, k], family=family)
    param.updated[, k] <- coefficients(fit)

    # model specific estimation of log evidence
    if (family$family == 'gaussian') {
      le.updated[, k] <- dnorm(y, predict(fit), sigma, log=TRUE)
    } else if (family$family == 'poisson') {
      le.updated[, k] <- dpois(y, predict(fit), log=TRUE)
    } else if (family$family == 'binomial') {
      le.updated[, k] <- dbinom(y, size = 1, prob = predict(fit), log = TRUE)
    } else if (family$family == 'Gamma') {
      le.updated[, k] = dgamma(y, shape = predict(fit), scale = predict(fit), log = TRUE)
    }
  }
  return(list(model = family, n = n, loglik = loglik, post = post, post.cp = post.cp, cp = cp, le.updated = le.updated, param.before = beta, param.updated = param.updated, response.variable = y))
}

logsumexp=function(l)
{
  i <- which.max(l);
  res <- l[i] + log1p(sum(exp(l[-i] - l[i])));
  if (is.nan(res)) res <- -Inf;
  return(res);
}
