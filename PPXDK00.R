# ===== momentum strategy on SPX daily data =====
# 1️⃣ Read CSV
file_path <- "/Users/joe/Downloads/spx.csv"
df <- read.csv(file_path, stringsAsFactors = FALSE)

# 2️⃣ Convert the Date column properly
df$Date <- as.Date(df$Date, format = "%m/%d/%Y")

# 3️⃣ Remove commas from numeric columns (they’re characters right now)
df$Open  <- as.numeric(gsub(",", "", df$Open))
df$High  <- as.numeric(gsub(",", "", df$High))
df$Low   <- as.numeric(gsub(",", "", df$Low))
df$Close <- as.numeric(gsub(",", "", df$Close))

# 4️⃣ Check types
str(df)
rf_annual <- 0.00
cost_per_turn <- 0.0005
ppy <- 252                   # daily data
lookback <- 20               # 20-day momentum window
top_k_drawdowns <- 5

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
  equity <- cumprod(1 + pnl)
  years <- as.numeric(max(dates) - min(dates)) / 365.25
  (tail(equity, 1))^(1/years) - 1
}
ann_vol <- function(pnl, ppy) sd(pnl, na.rm=TRUE) * sqrt(ppy)
sharpe <- function(pnl, ppy, rf_annual) {
  vol <- ann_vol(pnl, ppy)
  mu <- mean(pnl, na.rm=TRUE) * ppy
  if(vol > 0) (mu - rf_annual)/vol else NA_real_
}
max_dd_fn <- function(pnl) {
  eq <- cumprod(1+pnl); peaks <- cummax(eq); min(eq/peaks - 1, na.rm=TRUE)
}
win_rate <- function(pnl) mean(pnl>0, na.rm=TRUE)

# --- Load and clean data ---
file_path <- "/Users/joe/Downloads/spx.csv"
df <- read.csv(file_path, stringsAsFactors = FALSE)

# Convert date first
df$Date <- as.Date(df$Date, format = "%m/%d/%Y")

# Order and remove duplicates
df <- df[order(df$Date), ]
df <- df[!duplicated(df$Date), ]

# --- Force numeric columns to stay numeric ---
num_cols <- c("Open", "High", "Low", "Close")
df[num_cols] <- lapply(df[num_cols], function(x) {
  as.numeric(gsub(",", "", x))
})

# Verify
str(df)
head(df)

# --- Compute returns ---
df$logP <- log(df$Close)
ret <- c(NA, diff(df$logP))
ret_dates <- df$Date[-1]

# --- Parameters ---
lookback <- 20  # lookback window (e.g., 20 days or 1 month if daily data)

# --- Compute log prices ---
df$logP <- log(df$Close)

# --- Initialize vector ---
mom <- rep(NA_real_, nrow(df))

# --- Compute momentum safely ---
for (i in seq_along(df$logP)) {
  if (i - lookback - 1 > 0) {
    mom[i] <- df$logP[i - 1] - df$logP[i - 1 - lookback]
  } else {
    mom[i] <- NA_real_
  }
}

df$momentum <- mom
head(df, 25)

# --- Generate positions ---
pos <- rep(0, length(mom))
for(i in seq_along(mom)){
  if(is.na(mom[i])) pos[i] <- 0
  else if(mom[i] > 0) pos[i] <- 1
  else if(mom[i] < 0) pos[i] <- -1
}

# --- Strategy returns ---
pos_lag <- c(0, pos[-length(pos)])
turnover <- c(0, abs(diff(pos)))
costs <- turnover * cost_per_turn
strat_ret <- pos_lag * ret - costs
strat_ret <- strat_ret[-1]
strat_dates <- ret_dates[-1]

# --- Metrics ---
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

cat("========== SPX Momentum Strategy ==========\n")
cat(sprintf("Lookback=%d days, Cost=%.4f\n", lookback, cost_per_turn))
cat(sprintf("Annual Return (CAGR): %s\n", fmt_pct(str_ann_ret)))
cat(sprintf("Annual Volatility   : %s\n", fmt_pct(str_ann_vol)))
cat(sprintf("Sharpe Ratio        : %.3f\n", str_sharpe))
cat(sprintf("Max Drawdown        : %s\n", fmt_pct(str_mdd)))
cat(sprintf("Win Rate            : %s\n", fmt_pct(str_wr)))

