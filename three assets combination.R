# -------------------------- Config ------------------------------
files <- list(
  PCAAQ03 = "PCAAQ03.csv",
  PPXDM00 = "PPXDM00.csv",
  PPXDK00 = "PPXDK00.csv"
)

ppy              <- 12          # periods per year (monthly)
trade_cost       <- 0.0000      # per position change
use_winsor       <- TRUE        # winsorize single-period returns
winsor_p         <- 0.01        # 1% tails
spread_scale_win <- 12          # months for spread Δ scaling
clamp_spread     <- 0.30        # cap |spread return| per month
auto_flip_sign   <- TRUE        # flip signal if Sharpe<0 (mean-reversion)

# MA search grid (per asset)
grid_short <- c(2,3,4,5,6)
grid_long  <- c(8,10,12,15,18,21,24)

# ------------------------ Utilities -----------------------------
lag_vec   <- function(x,k=1) if (k<=0) x else c(rep(NA,k), head(x,-k))

roll_mean <- function(x, n){
  x <- as.numeric(x)
  if (is.na(n) || n <= 1) return(x)
  if (length(x) < n || sum(is.finite(x)) < n) return(rep(NA_real_, length(x)))
  as.numeric(stats::filter(x, rep(1/n, n), sides = 1))
}

fmt_pct <- function(x) ifelse(is.na(x), "NA", sprintf("%.2f%%", 100*x))

ann_vol <- function(r, ppy) stats::sd(r, na.rm=TRUE) * sqrt(ppy)

sharpe  <- function(r, ppy, rf=0){
  v  <- ann_vol(r, ppy)
  mu <- mean(r, na.rm=TRUE) * ppy
  if (!is.finite(v) || v == 0) return(NA_real_) else (mu - rf)/v
}

max_dd <- function(r){
  r <- r[is.finite(r)]
  if (length(r) < 2) return(0)
  eq <- suppressWarnings(cumprod(1 + r))
  eq <- eq[is.finite(eq) & eq > 0]
  if (length(eq) < 2) return(0)
  peaks <- cummax(eq)
  min(eq/peaks - 1, na.rm = TRUE)
}

cagr_from_r <- function(r, ppy){
  r <- r[is.finite(r)]
  if (length(r) < 2) return(NA_real_)
  eq  <- cumprod(1 + r)
  yrs <- length(r) / ppy
  if (yrs <= 0) return(NA_real_)
  tail(eq, 1)^(1/yrs) - 1
}

# ---------------------- CSV reader (robust) ---------------------
read_platts_csv <- function(path){
  con <- file(path, "r", encoding="UTF-8"); on.exit(close(con))
  first <- readLines(con, n=1, warn=FALSE)
  skip <- ifelse(tolower(trimws(first)) == "sep=,", 1, 0)
  df <- read.csv(path, stringsAsFactors=FALSE, check.names=FALSE, skip=skip)
  
  # detect date column
  date_col <- NULL
  for (nm in names(df)) if (nm %in% c("Timestamp","Date","date","DATE","timestamp")) { date_col <- nm; break }
  if (is.null(date_col)){
    for (nm in names(df)) {
      v <- suppressWarnings(as.Date(df[[nm]]))
      if (mean(!is.na(v)) > 0.8) { date_col <- nm; break }
    }
  }
  if (is.null(date_col)) stop("Date column not found in: ", path)
  
  # detect price/series column
  price_col <- NULL
  tag <- names(df)[grepl(": CLOSE", names(df), fixed=TRUE)]
  if (length(tag) > 0) price_col <- tag[1]
  if (is.null(price_col)){
    cand <- intersect(c("Close","close","Adj Close","Adj_Close","Price","price","VALUE"), names(df))
    if (length(cand) > 0) price_col <- cand[1]
  }
  if (is.null(price_col)){
    nums <- names(df)[sapply(df, function(x) suppressWarnings(!all(is.na(as.numeric(gsub(",","",as.character(x)))))))]
    if (length(nums) > 0) price_col <- tail(nums, 1)
  }
  if (is.null(price_col)) stop("Price column not found in: ", path)
  
  Date  <- suppressWarnings(as.Date(df[[date_col]]))
  Price <- suppressWarnings(as.numeric(gsub(",","",as.character(df[[price_col]]))))
  x <- data.frame(Date=Date, Price=Price, stringsAsFactors=FALSE)
  x <- x[!is.na(x$Date) & !is.na(x$Price), ]
  x <- x[order(x$Date), ]
  x <- x[!duplicated(x$Date), ]
  
  # end-of-month sample (last obs per YYYY-MM)
  if (nrow(x) > 0){
    ym <- format(x$Date, "%Y-%m")
    idx_last <- tapply(seq_len(nrow(x)), ym, tail, 1)
    x <- x[unlist(idx_last), ]
    rownames(x) <- NULL
  }
  x
}

# ------------------ Return definitions --------------------------
# price assets: log returns
ret_price <- function(p){
  p <- as.numeric(p); p[p <= 0] <- NA
  r <- c(NA, diff(log(p)))
  if (use_winsor){
    fr <- r[is.finite(r)]
    if (length(fr) >= 10){
      lo <- stats::quantile(fr, winsor_p, na.rm=TRUE)
      hi <- stats::quantile(fr, 1 - winsor_p, na.rm=TRUE)
      r  <- pmin(pmax(r, lo), hi)
    }
  }
  r
}

# spread assets: scaled differences (Δ/rolling scale) with clamp
ret_spread <- function(p){
  p <- as.numeric(p)
  L <- spread_scale_win
  scale <- roll_mean(abs(p), L)
  d  <- c(NA, diff(p))
  r  <- d / (scale + 1e-8)
  r  <- pmax(pmin(r, clamp_spread), -clamp_spread)
  if (use_winsor){
    fr <- r[is.finite(r)]
    if (length(fr) >= 10){
      lo <- stats::quantile(fr, winsor_p, na.rm=TRUE)
      hi <- stats::quantile(fr, 1 - winsor_p, na.rm=TRUE)
      r  <- pmin(pmax(r, lo), hi)
    }
  }
  r
}

# ------------------ Signal from MA cross ------------------------
make_signal <- function(price, short_n, long_n){
  sma_s <- roll_mean(price, short_n)
  sma_l <- roll_mean(price, long_n)
  s_now <- ifelse(sma_s > sma_l,  1L, ifelse(sma_s < sma_l, -1L, 0L))
  s_pos <- lag_vec(s_now, 1); s_pos[is.na(s_pos)] <- 0
  s_pos
}

# ------------------ Grid search per asset -----------------------
pick_best_ma <- function(df_price, df_ret, asset_name){
  best <- list(sh=-Inf, s=NA_integer_, l=NA_integer_, flip=FALSE)
  n_eff <- sum(is.finite(df_price))
  for (s in grid_short){
    for (l in grid_long){
      if (l <= s + 2) next
      if (l > n_eff) next
      sig <- make_signal(df_price, s, l)
      if (all(!is.finite(sig))) next
      
      pnl <- sig * df_ret
      sc  <- sharpe(pnl, ppy, 0)
      if (is.finite(sc) && sc > best$sh){
        best <- list(sh=sc, s=s, l=l, flip=FALSE)
      }
      pnl2 <- (-sig) * df_ret
      sc2  <- sharpe(pnl2, ppy, 0)
      if (is.finite(sc2) && sc2 > best$sh){
        best <- list(sh=sc2, s=s, l=l, flip=TRUE)
      }
    }
  }
  best
}

# ------------------ Prepare each asset --------------------------
raw <- lapply(files, read_platts_csv)

prep_per_asset <- list()
for (nm in names(raw)){
  df <- raw[[nm]]
  ret <- if (nm == "PPXDM00") ret_spread(df$Price) else ret_price(df$Price)
  
  best <- pick_best_ma(df$Price, ret, nm)
  sig  <- make_signal(df$Price, best$s, best$l)
  if (auto_flip_sign && isTRUE(best$flip)) sig <- -sig
  
  prep_per_asset[[nm]] <- data.frame(
    Date    = df$Date,
    YM      = format(df$Date, "%Y-%m"),
    ret     = ret,
    sig     = sig,
    s       = best$s,
    l       = best$l,
    flipped = best$flip,
    stringsAsFactors = FALSE
  )
}

# ------------------ Align by YearMonth --------------------------
merge_two <- function(a,b){
  merge(a[,c("YM","ret","sig")], b[,c("YM","ret","sig")],
        by="YM", all=FALSE, suffixes=c(".x",".y"))
}
tmp    <- merge_two(prep_per_asset[[1]], prep_per_asset[[2]])
merged <- merge(tmp, prep_per_asset[[3]][,c("YM","ret","sig")], by="YM", all=FALSE)
names(merged) <- c("YM",
                   paste0(names(prep_per_asset)[1], c("_ret","_sig")),
                   paste0(names(prep_per_asset)[2], c("_ret","_sig")),
                   paste0(names(prep_per_asset)[3], c("_ret","_sig")))
merged <- merged[complete.cases(merged), ]
stopifnot(nrow(merged) >= 12)

# ------------------ Single-asset PnL ----------------------------
single_pnls <- list()
for (nm in names(prep_per_asset)){
  r   <- merged[[paste0(nm,"_ret")]]
  sig <- merged[[paste0(nm,"_sig")]]
  pos <- sig
  turn <- c(0, abs(diff(pos)))
  pnl <- pos * r - trade_cost * turn
  single_pnls[[nm]] <- pnl
}

# ------------------ Portfolio constructions ---------------------
ret_mat <- as.matrix(merged[, paste0(names(prep_per_asset), "_ret"), drop=FALSE])
ew_ret  <- rowMeans(ret_mat, na.rm = FALSE)

sig_mat <- as.matrix(merged[, paste0(names(prep_per_asset), "_sig"), drop=FALSE])
sig_sum <- rowSums(sig_mat, na.rm=TRUE)

# Consensus on equal-weight returns
cons_pos  <- ifelse(sig_sum >= 2,  1L, ifelse(sig_sum <= -2, -1L, 0L))
cons_turn <- c(0, abs(diff(cons_pos)))
cons_pnl  <- cons_pos * ew_ret - trade_cost * cons_turn

# Weighted-signal portfolio (inverse-vol weights)
vols <- apply(ret_mat, 2, function(x) stats::sd(x, na.rm=TRUE))
w_iv <- 1/vols; w_iv[!is.finite(w_iv)] <- 0
w_iv <- w_iv / sum(w_iv)
sig_w  <- sweep(sig_mat, 2, w_iv, `*`)
ws_turn <- c(0, rowSums(abs(diff(sig_w)), na.rm=TRUE))
ws_pnl  <- rowSums(sig_w * ret_mat, na.rm=TRUE) - trade_cost * ws_turn

# ------------------ Metrics & Summary ---------------------------
mk_metrics <- function(r){
  list(
    Annual_Return = cagr_from_r(r, ppy),
    Annual_Vol    = ann_vol(r, ppy),
    Sharpe        = sharpe(r, ppy, 0),
    Max_Drawdown  = max_dd(r),
    Win_Rate      = mean(r > 0, na.rm=TRUE)
  )
}

res <- list(
  PCAAQ03 = mk_metrics(single_pnls$PCAAQ03),
  PPXDM00 = mk_metrics(single_pnls$PPXDM00),
  PPXDK00 = mk_metrics(single_pnls$PPXDK00),
  Consensus_EW_Portfolio = mk_metrics(cons_pnl),
  WeightedSignal_IV      = mk_metrics(ws_pnl)
)

summ <- data.frame(
  Strategy       = names(res),
  Annual_Return  = sapply(res, `[[`, "Annual_Return"),
  Annual_Vol     = sapply(res, `[[`, "Annual_Vol"),
  Sharpe         = sapply(res, `[[`, "Sharpe"),
  Max_Drawdown   = sapply(res, `[[`, "Max_Drawdown"),
  Win_Rate       = sapply(res, `[[`, "Win_Rate"),
  row.names = NULL
)

summ_fmt <- summ
for (cc in c("Annual_Return","Annual_Vol","Max_Drawdown","Win_Rate"))
  summ_fmt[[cc]] <- sapply(summ_fmt[[cc]], fmt_pct)

cat("========== Three-Asset Momentum (Auto-Selected MA, robust returns) ==========\n")
cat(sprintf("cost=%.4f, winsor=%s @%.1f%%, spread_scale_win=%d, clamp=±%.0f%%, auto_flip=%s\n\n",
            trade_cost, use_winsor, 100*winsor_p, spread_scale_win, 100*clamp_spread, auto_flip_sign))
for (nm in names(prep_per_asset)){
  s <- unique(prep_per_asset[[nm]]$s); l <- unique(prep_per_asset[[nm]]$l); f <- unique(prep_per_asset[[nm]]$flipped)
  cat(sprintf("Chosen MAs for %-8s : short=%d, long=%d, flipped=%s\n", nm, s, l, f))
}
cat("\n")
print(summ_fmt, row.names = FALSE)
