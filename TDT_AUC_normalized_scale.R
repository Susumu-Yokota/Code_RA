# -----------------------------
# 割引課題 正規化AUC計算
# -----------------------------

# パッケージ（必要なら）
# install.packages("pracma")
library(pracma)

# データ読み込み
data <- read.csv("discounting.csv")

# 列名を扱いやすく（必要に応じて調整）
# 例: subject, delay, amount, V
colnames(data) <- c("ID", "delay", "amount", "V")

# AUC計算（正規化込み）
auc_results <- data.frame(
  ID = unique(data$ID),
  AUC = sapply(unique(data$ID), function(id) {
    
    df <- data[data$ID == id, ]
    
    # ---- 並び替え（超重要）----
    df <- df[order(df$delay), ]
    
    # ---- 正規化 ----
    x <- df$delay / max(df$delay)   # delayを0-1に
    y <- df$V / df$amount           # valueを0-1に（今回は1000で割るのと同じ）
    
    # ---- AUC計算 ----
    pracma::trapz(x, y)
  })
)

# z変換（そのままでもOK）
auc_results$AUC_z <- as.numeric(scale(auc_results$AUC))

# 保存
write.csv(auc_results, "AUC_results_normalized.csv", row.names = FALSE)