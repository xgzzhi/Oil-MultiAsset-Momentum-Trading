# ===== minimal mean-reversion pipeline for your CSV shape =====

file_path <- "PCAAQ03.csv"
rf_annual <- 0.00
roll_window <- 60          # rolling window for z-score (periods = months here)
z_entry <- 1.0
z_exit  <- 0.2
cost_per_turn <- 0.0005
top_k_drawdowns <- 5
ppy <- 12                  # monthly data

fmt_pct <- function(x) sprintf("%.2f%%", 100*x)

compute_drawdowns <- function(equity, dates){
  peaks <- cummax(equity)
  res <- data.frame(PeakDate=as.Date(character()),TroughDate=as.Date(character()),
                    RecoveryDate=as.Date(character()),Depth=numeric(),
                    Duration=integer(),ToTrough=integer(),ToRecovery=integer(),
                    stringsAsFactors=FALSE)
  i<-1; n<-length(equity)
  while(i<n){
    if(equity[i]<peaks[i]){
      peak_val <- peaks[i]
      peak_pos <- max(which(peaks[1:i]==peak_val))
      j<-i; trough_pos<-i
      while(j<=n && equity[j]<peak_val){
        if(equity[j]<=equity[trough_pos]) trough_pos<-j
        j<-j+1
      }
      recovery_pos <- if(j<=n) j else NA_integer_
      depth <- equity[trough_pos]/peak_val - 1
      dur <- if(is.na(recovery_pos)) (n-peak_pos) else (recovery_pos-peak_pos)
      to_trough <- trough_pos-peak_pos
      to_recovery <- if(is.na(recovery_pos)) NA_integer_ else (recovery_pos-trough_pos)
      res <- rbind(res, data.frame(
        PeakDate=dates[peak_pos], TroughDate=dates[trough_pos],
        RecoveryDate=if(is.na(recovery_pos)) as.Date(NA) else dates[recovery_pos],
        Depth=depth, Duration=dur, ToTrough=to_trough, ToRecovery=to_recovery,
        stringsAsFactors=FALSE))
      i <- if(is.na(recovery_pos)) n else recovery_pos
    } else i <- i+1
  }
  res[order(res$Depth), , drop=FALSE]
}

cagr <- function(pnl, dates) {
  if(length(pnl)==0) return(NA_real_)
  equity <- cumprod(1 + pnl)
  years <- as.numeric(max(dates) - min(dates))/365.25
  if(years<=0) return(NA_real_)
  (tail(equity,1))^(1/years) - 1
}
ann_vol <- function(pnl, ppy) stats::sd(pnl, na.rm=TRUE)*sqrt(ppy)
sharpe <- function(pnl, ppy, rf_annual){
  vol <- ann_vol(pnl, ppy); mu <- mean(pnl, na.rm=TRUE)*ppy
  if(vol>0) (mu - rf_annual)/vol else NA_real_
}
max_dd_fn <- function(pnl){
  eq <- cumprod(1+pnl); peaks <- cummax(eq); min(eq/peaks - 1, na.rm=TRUE)
}
win_rate <- function(pnl) mean(pnl>0, na.rm=TRUE)

# --- read your file (skip the "sep=," hint line) ---
df_raw <- read.csv(file_path, stringsAsFactors=FALSE, check.names=FALSE, skip=1)

# pick columns
date_col  <- "Timestamp"
price_col <- grep(": CLOSE", names(df_raw), value=TRUE)[1]

# make Date & Price
df <- data.frame(
  Date  = as.Date(df_raw[[date_col]], format="%Y-%m-%d"),
  Price = as.numeric(gsub(",", "", df_raw[[price_col]])),
  stringsAsFactors = FALSE
)
df <- df[!is.na(df$Date) & !is.na(df$Price), ]
df <- df[order(df$Date), ]
df <- df[!duplicated(df$Date), ]
stopifnot(nrow(df) >= 3)

# --- series ---
logP <- log(df$Price)
ret  <- c(NA, diff(logP))
ret   <- ret[!is.na(ret)]
ret_dates <- df$Date[-1]
years_span <- as.numeric(max(df$Date) - min(df$Date))/365.25

# --- OU/AR(1) estimation on logP increments ---
Xlag <- logP[-length(logP)]
dX   <- diff(logP)
ols  <- lm(dX ~ Xlag)
a <- coef(ols)[1]; b <- coef(ols)[2]
kappa <- -as.numeric(b)
mu    <- -as.numeric(a/b)
halflife_periods <- if (kappa>0) log(2)/kappa else NA_real_
halflife_years   <- if (is.na(halflife_periods)) NA_real_ else halflife_periods/ppy

cat("========== OU/AR(1) Estimation ==========\n")
cat(sprintf("a = %.6f, b = %.6f  =>  kappa = %.6f,  mu(logP) = %.6f\n", a, b, kappa, mu))
cat(sprintf("Half-life: %.2f periods  (~ %.2f years @ ppy=%d)\n\n",
            halflife_periods, halflife_years, ppy))

# --- z-score mean-reversion strategy on logP ---
roll_mean <- rep(NA_real_, length(logP))
roll_sd   <- rep(NA_real_, length(logP))
for(i in seq_along(logP)){
  j0 <- max(1, i-roll_window+1)
  if(i-j0+1 >= 5){
    roll_mean[i] <- mean(logP[j0:i], na.rm=TRUE)
    roll_sd[i]   <- stats::sd(logP[j0:i], na.rm=TRUE)
  }
}
z <- (logP - roll_mean)/roll_sd

pos <- rep(0, length(z))
for(i in seq_along(z)){
  if(is.na(z[i])) { pos[i] <- if(i>1) pos[i-1] else 0; next }
  if(abs(z[i]) < z_exit) {
    pos[i] <- 0
  } else if(z[i] <= -z_entry){
    pos[i] <- 1
  } else if(z[i] >=  z_entry){
    pos[i] <- -1
  } else {
    pos[i] <- if(i>1) pos[i-1] else 0
  }
}
turnover <- c(0, abs(diff(pos)))
costs <- turnover * cost_per_turn
pos_lag <- c(0, pos[-length(pos)])
strat_ret <- pos_lag * ret - costs[-1]
strat_dates <- ret_dates

# --- metrics: underlying vs strategy ---
under_ann_ret <- cagr(ret[-1], ret_dates[-1])
under_ann_vol <- ann_vol(ret, ppy)
under_sharpe  <- sharpe(ret, ppy, rf_annual)
under_mdd     <- max_dd_fn(ret)
under_wr      <- win_rate(ret)

str_ann_ret <- cagr(strat_ret, strat_dates)
str_ann_vol <- ann_vol(strat_ret, ppy)
str_sharpe  <- sharpe(strat_ret, ppy, rf_annual)
str_mdd     <- max_dd_fn(strat_ret)
str_wr      <- win_rate(strat_ret)

equity_str <- cumprod(1 + strat_ret)
dd_tbl <- compute_drawdowns(c(1, equity_str), c(df$Date[1], strat_dates))
dd_top <- if(nrow(dd_tbl)>0) head(dd_tbl[order(dd_tbl$Depth),], top_k_drawdowns) else dd_tbl

cat("========== Frequency & Sample ==========\n")
cat(sprintf("Data Range : %s ~ %s (%.2f years, ppy=%d)\n\n",
            format(min(df$Date)), format(max(df$Date)), years_span, ppy))

cat("========== Underlying (Buy&Hold) ==========\n")
cat(sprintf("Annual Return (CAGR): %s\n", fmt_pct(under_ann_ret)))
cat(sprintf("Annual Volatility   : %s\n", fmt_pct(under_ann_vol)))
cat(sprintf("Sharpe (rf=%.2f%%)  : %.3f\n", 100*rf_annual, under_sharpe))
cat(sprintf("Max Drawdown        : %s\n", fmt_pct(under_mdd)))
cat(sprintf("Win Rate            : %s\n\n", fmt_pct(under_wr)))

cat("========== Mean-Reversion Strategy ==========\n")
cat(sprintf("Params: window=%d, z_entry=%.2f, z_exit=%.2f, cost=%.4f\n",
            roll_window, z_entry, z_exit, cost_per_turn))
cat(sprintf("Annual Return (CAGR): %s\n", fmt_pct(str_ann_ret)))
cat(sprintf("Annual Volatility   : %s\n", fmt_pct(str_ann_vol)))
cat(sprintf("Sharpe (rf=%.2f%%)  : %.3f\n", 100*rf_annual, str_sharpe))
cat(sprintf("Max Drawdown        : %s\n", fmt_pct(str_mdd)))
cat(sprintf("Win Rate            : %s\n\n", fmt_pct(str_wr)))

cat("========== Strategy Top Drawdown Patterns ==========\n")
if(nrow(dd_top)==0){
  cat("No drawdown segments (or insufficient data).\n")
} else {
  dd_out <- dd_top; dd_out$Depth <- fmt_pct(dd_out$Depth)
  print(dd_out, row.names=FALSE)
}


