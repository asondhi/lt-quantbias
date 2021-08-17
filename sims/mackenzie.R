# Code from Mackenzie paper

dep.truncation.Survival <- function(time, truncation.time, status, Prob.Trunc=NULL,
                                    Lower.Bound=NULL) {
  keep <- !is.na(time) & !is.na(status) & !is.na(truncation.time)
  t <- time[keep]
  v <- truncation.time[keep]
  s <- status[keep]
  if (is.null(Lower.Bound)) Lower.Bound <- 0
  ord <- order(v)
  v <- v[ord]
  t <- t[ord]
  s <- s[ord]
  o.cox <- coxph(Surv(v, t, s) ~ v)
  lin <- as.matrix(v) %*% o.cox$coef
  lin <- lin - mean(lin)
  o.CH <- Cumulative.Hazard(t, v, s, lin)
  n.t <- length(o.CH$time)
  tt <- c(0, o.CH$time)
  CH <- c(0, o.CH$Cum.Haz)
  if (is.null(Prob.Trunc)) {
    Prob.Trunc <- rep(NA, length(v))
    i.t <- 1+n.t
    for (i in length(v):1) {
      while(tt[i.t] > v[i]) i.t <- i.t-1
      Prob.Trunc[i] <- exp(-exp(lin[i])*CH[i.t])
    }
  }
  
  if (!is.null(Prob.Trunc)) {
    Prob.Trunc <- (Prob.Trunc[keep])[ord]
  }
  # Use lower bound
  Prob.Trunc <- ifelse(Prob.Trunc<Lower.Bound, Lower.Bound, Prob.Trunc)
  Q <- 1/mean(1/Prob.Trunc)
  S <- rep(NA, n.t)
  for (i.t in 1:n.t) {
    S[i.t] <- Q * mean(exp(-exp(lin)*o.CH$Cum.Haz[i.t]) / Prob.Trunc)
  }
  S.cond <- list(time=o.CH$time, surv=exp(-o.CH$Cum.Haz))
  list(survival=S, time=o.CH$time, Q=Q, log.HR=o.cox$coef, cox.iter=o.cox$iter,
       S.cond=S.cond)
}

Cumulative.Hazard <- function(time, time.start=0, status, x, correction.power=1/3)
{
  if (length(time.start)==0) time.start <- rep(0, length(time))
  n <- length(time)
  ot <- order(time)
  ti.ot <- time[ot]
  ti.start.ot <- time.start[ot]
  st.ot <- status[ot]
  x.ot <- x[ot]
  uniq.ev.ti <- unique(ti.ot[st.ot==1])
  n.uniq.ev.ti <- length(uniq.ev.ti)
  H <- rep(NA, n.uniq.ev.ti)
  nr <- rep(NA, n)
  for (i in 1:n.uniq.ev.ti) {
    Y.i.t <- (ti.ot>=uniq.ev.ti[i]) * (ti.start.ot < uniq.ev.ti[i])
    H[i] <- ifelse(sum(Y.i.t)<= n^correction.power, 0, sum(st.ot==1 & ti.ot==uniq.ev.ti[i])
                   /sum(Y.i.t * exp(x.ot)))
  }
  list(time=uniq.ev.ti, Cum.Haz = cumsum(H))
}

# Estimate proportion of LT, accounting for censoring
EstimatePropLT <- function(data_lt) {
  
  data_lt <- data_lt %>%
    mutate(ones = 1)
  
  m_entry <- coxph(Surv(entry, time, ones) ~ entry, data = data_lt)
  bh <- basehaz(m_entry, centered = FALSE)
  bh_fn <- stepfun(bh$time, c(0, bh$hazard))
  surv_prob <- exp(-bh_fn(data_lt$entry) * exp(predict(m_entry, type = "lp")))
  
  Q_est <- 1 / mean(1 / surv_prob)
  
  return(1 - Q_est)
}

# Estimate median survival
EstimateMedianSurvival <- function(data_lt) {
  
  m_entry <- coxph(Surv(entry, time, status) ~ entry, data = data_lt)
  
  mkdep <- dep.truncation.Survival(time = data_lt$time, 
                                   truncation.time = data_lt$entry, 
                                   status = data_lt$status)
  times_obs <- mkdep$time
  surv_prob_est <- mkdep$survival

  median_est <- times_obs[which.min(abs(surv_prob_est - 0.5))]

  return(median_est)
}


  