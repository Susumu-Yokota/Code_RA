# ============================================================
# 階層的重回帰（Composite + Interaction）
# ASD特性中心（AQ）
# Step4a: AQ×EF
# Step4b: AQ×DM
# Step4c: EF×DM
# 全StepのR2 / ΔR2 / F検定 + 最終モデルβ
# 有意 outcome の simple slope
# ============================================================

library(car)
library(lm.beta)
library(ggplot2)

# ---------- データ ----------
df <- read.csv("regression_composit.csv")

# ---------- 目的変数 ----------
#outcomes <- c(
#  "Fri_Com_A_fair","Fri_Com_B_fair","Fri_Com_C_fair","Fri_Com_D_fair",
#  "Fri_Oth_A_fair","Fri_Oth_B_fair","Fri_Oth_C_fair","Fri_Oth_D_fair",
#  "Oth_Com_A_fair","Oth_Com_B_fair","Oth_Com_C_fair","Oth_Com_D_fair"
#)

outcomes <- c(
  "ComA_fair","ComB_fair","ComC_fair","ComD_fair",
  "FriA_fair","FriB_fair","FriC_fair","FriD_fair",
  "OthA_fair","OthB_fair","OthC_fair","OthD_fair"
)

# ============================================================
# 回帰
# ============================================================

results <- lapply(outcomes, function(y) {
  
  # ----- モデル -----
  m0 <- lm(as.formula(paste(y, "~ age + sex")), data = df)
  m1 <- lm(as.formula(paste(y, "~ age + sex + AQ_total_Z")), data = df)
  m2 <- lm(as.formula(paste(y, "~ age + sex + AQ_total_Z + EF_composit_Z")), data = df)
  m3 <- lm(as.formula(paste(y, "~ age + sex + AQ_total_Z + EF_composit_Z + DM_composit_Z")), data = df)
  
  m4a <- lm(as.formula(paste(
    y, "~ age + sex + AQ_total_Z + EF_composit_Z + DM_composit_Z + AQ_total_Z:EF_composit_Z"
  )), data = df)
  
  m4b <- lm(as.formula(paste(
    y, "~ age + sex + AQ_total_Z + EF_composit_Z + DM_composit_Z + AQ_total_Z:DM_composit_Z"
  )), data = df)
  
  m4c <- lm(as.formula(paste(
    y, "~ age + sex + AQ_total_Z + EF_composit_Z + DM_composit_Z + EF_composit_Z:DM_composit_Z"
  )), data = df)
  
  # ----- R2 -----
  R2_0 <- summary(m0)$r.squared
  R2_1 <- summary(m1)$r.squared
  R2_2 <- summary(m2)$r.squared
  R2_3 <- summary(m3)$r.squared
  R2_4a <- summary(m4a)$r.squared
  R2_4b <- summary(m4b)$r.squared
  R2_4c <- summary(m4c)$r.squared
  
  # ----- ΔR2 -----
  d01 <- anova(m0, m1)
  d12 <- anova(m1, m2)
  d23 <- anova(m2, m3)
  d34a <- anova(m3, m4a)
  d34b <- anova(m3, m4b)
  d34c <- anova(m3, m4c)
  
  # ----- F検定 -----
  f <- function(m) summary(m)$fstatistic
  f0 <- f(m0); f1 <- f(m1); f2 <- f(m2); f3 <- f(m3)
  f4a <- f(m4a); f4b <- f(m4b); f4c <- f(m4c)
  
  # ----- VIF（多重共線性チェック） -----
  vif_m3 <- vif(m3)
  vif_m4a <- vif(m4a)
  vif_m4b <- vif(m4b)
  vif_m4c <- vif(m4c)
  
  # 最大VIFだけ保存（シンプル版）
  max_vif_m3  <- max(vif_m3)
  max_vif_m4a <- max(vif_m4a)
  max_vif_m4b <- max(vif_m4b)
  max_vif_m4c <- max(vif_m4c)
  
  # ----- β（Step4） -----
  b4a <- summary(lm.beta(m4a))$coefficients
  b4b <- summary(lm.beta(m4b))$coefficients
  b4c <- summary(lm.beta(m4c))$coefficients
  
  data.frame(
    outcome = y,
    
    # ---- R2 ----
    R2_0, R2_1, R2_2, R2_3, R2_4a, R2_4b, R2_4c,
    
    # ---- VIF ----
    max_vif_m3,
    max_vif_m4a,
    max_vif_m4b,
    max_vif_m4c,
    
    # ---- ΔR2 ----
    dR2_01 = R2_1 - R2_0, p_dR2_01 = d01$`Pr(>F)`[2],
    dR2_12 = R2_2 - R2_1, p_dR2_12 = d12$`Pr(>F)`[2],
    dR2_23 = R2_3 - R2_2, p_dR2_23 = d23$`Pr(>F)`[2],
    dR2_4a = R2_4a - R2_3, p_dR2_4a = d34a$`Pr(>F)`[2],
    dR2_4b = R2_4b - R2_3, p_dR2_4b = d34b$`Pr(>F)`[2],
    dR2_4c = R2_4c - R2_3, p_dR2_4c = d34c$`Pr(>F)`[2],
    
    # ---- F検定 ----
    F_0    = unname(f0["value"]),
    df1_0  = unname(f0["numdf"]),
    df2_0  = unname(f0["dendf"]),
    p_F_0  = pf(f0["value"], f0["numdf"], f0["dendf"], lower.tail = FALSE),
    
    F_1    = unname(f1["value"]),
    df1_1  = unname(f1["numdf"]),
    df2_1  = unname(f1["dendf"]),
    p_F_1  = pf(f1["value"], f1["numdf"], f1["dendf"], lower.tail = FALSE),
    
    F_2    = unname(f2["value"]),
    df1_2  = unname(f2["numdf"]),
    df2_2  = unname(f2["dendf"]),
    p_F_2  = pf(f2["value"], f2["numdf"], f2["dendf"], lower.tail = FALSE),
    
    F_3    = unname(f3["value"]),
    df1_3  = unname(f3["numdf"]),
    df2_3  = unname(f3["dendf"]),
    p_F_3  = pf(f3["value"], f3["numdf"], f3["dendf"], lower.tail = FALSE),
    
    F_4a   = unname(f4a["value"]),
    df1_4a = unname(f4a["numdf"]),
    df2_4a = unname(f4a["dendf"]),
    p_F_4a = pf(f4a["value"], f4a["numdf"], f4a["dendf"], lower.tail = FALSE),
    
    F_4b   = unname(f4b["value"]),
    df1_4b = unname(f4b["numdf"]),
    df2_4b = unname(f4b["dendf"]),
    p_F_4b = pf(f4b["value"], f4b["numdf"], f4b["dendf"], lower.tail = FALSE),
    
    F_4c   = unname(f4c["value"]),
    df1_4c = unname(f4c["numdf"]),
    df2_4c = unname(f4c["dendf"]),
    p_F_4c = pf(f4c["value"], f4c["numdf"], f4c["dendf"], lower.tail = FALSE),
    
    
    # ---- β & p（Step4）----
    # ---- Step4a ----
    beta_AQ_4a = b4a["AQ_total_Z", "Standardized"],
    p_AQ_4a    = b4a["AQ_total_Z", "Pr(>|t|)"],
    
    beta_EF_4a = b4a["EF_composit_Z", "Standardized"],
    p_EF_4a    = b4a["EF_composit_Z", "Pr(>|t|)"],
    
    beta_DM_4a = b4a["DM_composit_Z", "Standardized"],
    p_DM_4a    = b4a["DM_composit_Z", "Pr(>|t|)"],
    
    beta_AQ_EF = b4a["AQ_total_Z:EF_composit_Z", "Standardized"],
    p_AQ_EF    = b4a["AQ_total_Z:EF_composit_Z", "Pr(>|t|)"],
    
    # ---- Step4b ----
    beta_AQ_4b = b4b["AQ_total_Z", "Standardized"],
    p_AQ_4b    = b4b["AQ_total_Z", "Pr(>|t|)"],
    
    beta_EF_4b = b4b["EF_composit_Z", "Standardized"],
    p_EF_4b    = b4b["EF_composit_Z", "Pr(>|t|)"],
    
    beta_DM_4b = b4b["DM_composit_Z", "Standardized"],
    p_DM_4b    = b4b["DM_composit_Z", "Pr(>|t|)"],
    
    beta_AQ_DM = b4b["AQ_total_Z:DM_composit_Z", "Standardized"],
    p_AQ_DM    = b4b["AQ_total_Z:DM_composit_Z", "Pr(>|t|)"],
    
    # ---- Step4c ----
    beta_AQ_4c = b4c["AQ_total_Z", "Standardized"],
    p_AQ_4c    = b4c["AQ_total_Z", "Pr(>|t|)"],
    
    beta_EF_4c = b4c["EF_composit_Z", "Standardized"],
    p_EF_4c    = b4c["EF_composit_Z", "Pr(>|t|)"],
    
    beta_DM_4c = b4c["DM_composit_Z", "Standardized"],
    p_DM_4c    = b4c["DM_composit_Z", "Pr(>|t|)"],
    
    beta_EF_DM = b4c["EF_composit_Z:DM_composit_Z", "Standardized"],
    p_EF_DM    = b4c["EF_composit_Z:DM_composit_Z", "Pr(>|t|)"]
    
  )
})

results_df <- do.call(rbind, results)

# ---------- FDR補正（interactionのみ） ----------
results_df$p_AQ_EF_FDR <- p.adjust(results_df$p_AQ_EF, method = "fdr")
results_df$p_AQ_DM_FDR <- p.adjust(results_df$p_AQ_DM, method = "fdr")
results_df$p_EF_DM_FDR <- p.adjust(results_df$p_EF_DM, method = "fdr")

write.csv(
  results_df,
  "hierarchical_regression_AQ_interaction_full_results.csv",
  row.names = FALSE
)

# ============================================================
# 有意 outcome 抽出
# ============================================================

sig_AQ_EF <- results_df$outcome[
  results_df$p_dR2_4a < .05 & results_df$p_AQ_EF < .05
]

sig_AQ_DM <- results_df$outcome[
  results_df$p_dR2_4b < .05 & results_df$p_AQ_DM < .05
]

sig_EF_DM <- results_df$outcome[
  results_df$p_dR2_4c < .05 & results_df$p_EF_DM < .05
]

# ============================================================
# Simple slope（Step4モデル完全一致）
# ============================================================

plot_simple_slope <- function(y, form, focal, moderator, label_mod, fname) {
  
  model <- lm(form, data = df)
  
  focal_seq <- seq(min(df[[focal]]), max(df[[focal]]), length.out = 100)
  
  base <- data.frame(
    age = mean(df$age),
    sex = mean(df$sex),
    AQ_total_Z = mean(df$AQ_total_Z),
    EF_composit_Z = mean(df$EF_composit_Z),
    DM_composit_Z = mean(df$DM_composit_Z)
  )
  
  newdat <- do.call(rbind, lapply(c(-1,1), function(m){
    tmp <- base[rep(1,length(focal_seq)),]
    tmp[[focal]] <- focal_seq
    tmp[[moderator]] <- m
    tmp
  }))
  
  newdat$pred <- predict(model, newdata = newdat)
  newdat$mod_level <- factor(newdat[[moderator]],
                             labels = paste0(label_mod,c(" (-1SD)"," (+1SD)")))
  
  p <- ggplot(newdat,
              aes(x=.data[[focal]], y=pred, color=mod_level)) +
    geom_line(linewidth=1.2) +
    theme_minimal(base_size=14)
  
  ggsave(fname, p, width=6, height=4)
}

# ---------- AQ × EF ----------
for(y in sig_AQ_EF){
  plot_simple_slope(
    y,
    form = as.formula(paste(
      y,"~ age + sex + AQ_total_Z + EF_composit_Z + DM_composit_Z +
          AQ_total_Z:EF_composit_Z")),
    focal="EF_composit_Z",
    moderator="AQ_total_Z",
    label_mod="AQ",
    fname=paste0("simple_slope_AQ_EF_",y,".png")
  )
}

# ---------- AQ × DM ----------
for(y in sig_AQ_DM){
  plot_simple_slope(
    y,
    form = as.formula(paste(
      y,"~ age + sex + AQ_total_Z + EF_composit_Z + DM_composit_Z +
          AQ_total_Z:DM_composit_Z")),
    focal="DM_composit_Z",
    moderator="AQ_total_Z",
    label_mod="AQ",
    fname=paste0("simple_slope_AQ_DM_",y,".png")
  )
}

# ---------- EF × DM ----------
for(y in sig_EF_DM){
  plot_simple_slope(
    y,
    form = as.formula(paste(
      y,"~ age + sex + AQ_total_Z + EF_composit_Z + DM_composit_Z +
          EF_composit_Z:DM_composit_Z")),
    focal="EF_composit_Z",
    moderator="DM_composit_Z",
    label_mod="DM",
    fname=paste0("simple_slope_EF_DM_",y,".png")
  )
}

#============================================================
#  Simple slopeの統計量を出力
#============================================================

library(broom)

get_simple_slope <- function(model, focal, moderator, mod_value){
  
  coefs <- coef(model)
  vcov_mat <- vcov(model)
  
  # interaction項の名前（順序対応）
  int_name1 <- paste0(focal, ":", moderator)
  int_name2 <- paste0(moderator, ":", focal)
  
  int_name <- ifelse(int_name1 %in% names(coefs), int_name1, int_name2)
  
  # 傾き
  b <- coefs[focal] + coefs[int_name] * mod_value
  
  # SE計算
  var_b <- vcov_mat[focal, focal] +
    (mod_value^2) * vcov_mat[int_name, int_name] +
    2 * mod_value * vcov_mat[focal, int_name]
  
  se <- sqrt(var_b)
  
  t <- b / se
  df <- df.residual(model)
  p <- 2 * pt(-abs(t), df)
  
  return(data.frame(
    slope = b,
    SE = se,
    t = t,
    p = p,
    moderator_value = mod_value
  ))
}

#============================================================
#各interactionごとにsimple slopeを取得
#============================================================

simple_results <- list()

for(y in outcomes){
  
  # ---- AQ × EF ----
  model_4a <- lm(as.formula(paste(
    y,"~ age + sex + AQ_total_Z + EF_composit_Z + DM_composit_Z +
        AQ_total_Z:EF_composit_Z")), data=df)
  
  s1 <- get_simple_slope(model_4a, "EF_composit_Z", "AQ_total_Z", -1)
  s2 <- get_simple_slope(model_4a, "EF_composit_Z", "AQ_total_Z", 1)
  
  tmp <- rbind(s1, s2)
  tmp$outcome <- y
  tmp$interaction <- "AQ×EF"
  
  simple_results[[length(simple_results)+1]] <- tmp
  
  
  # ---- AQ × DM ----
  model_4b <- lm(as.formula(paste(
    y,"~ age + sex + AQ_total_Z + EF_composit_Z + DM_composit_Z +
        AQ_total_Z:DM_composit_Z")), data=df)
  
  s1 <- get_simple_slope(model_4b, "DM_composit_Z", "AQ_total_Z", -1)
  s2 <- get_simple_slope(model_4b, "DM_composit_Z", "AQ_total_Z", 1)
  
  tmp <- rbind(s1, s2)
  tmp$outcome <- y
  tmp$interaction <- "AQ×DM"
  
  simple_results[[length(simple_results)+1]] <- tmp
  
  
  # ---- EF × DM ----
  model_4c <- lm(as.formula(paste(
    y,"~ age + sex + AQ_total_Z + EF_composit_Z + DM_composit_Z +
        EF_composit_Z:DM_composit_Z")), data=df)
  
  s1 <- get_simple_slope(model_4c, "EF_composit_Z", "DM_composit_Z", -1)
  s2 <- get_simple_slope(model_4c, "EF_composit_Z", "DM_composit_Z", 1)
  
  tmp <- rbind(s1, s2)
  tmp$outcome <- y
  tmp$interaction <- "EF×DM"
  
  simple_results[[length(simple_results)+1]] <- tmp
}

simple_df <- do.call(rbind, simple_results)

write.csv(simple_df,
          "simple_slope_results.csv",
          row.names = FALSE)

#============================================================
#simple slopeで有意な傾きが得られなかった場合に、傾きの差の検定をする
#slope差（+1SD vs -1SD）の検定
#============================================================

get_slope_difference <- function(model, focal, moderator){
  
  coefs <- coef(model)
  vcov_mat <- vcov(model)
  
  # interaction項名
  int_name1 <- paste0(focal, ":", moderator)
  int_name2 <- paste0(moderator, ":", focal)
  int_name <- ifelse(int_name1 %in% names(coefs), int_name1, int_name2)
  
  # 傾き差：(+1SD) - (-1SD) = 2 * b_interaction
  diff_b <- 2 * coefs[int_name]
  
  # SE（分散共分散から）
  var_diff <- 4 * vcov_mat[int_name, int_name]
  se_diff <- sqrt(var_diff)
  
  t <- diff_b / se_diff
  df <- df.residual(model)
  p <- 2 * pt(-abs(t), df)
  
  return(data.frame(
    slope_diff = diff_b,
    SE = se_diff,
    t = t,
    p = p
  ))
}

#============================================================
#各interactionごとにslope差を取得
#============================================================
slope_diff_results <- list()

for(y in outcomes){
  
  # ---- AQ × EF ----
  model_4a <- lm(as.formula(paste(
    y,"~ age + sex + AQ_total_Z + EF_composit_Z + DM_composit_Z +
        AQ_total_Z:EF_composit_Z")), data=df)
  
  d <- get_slope_difference(model_4a, "EF_composit_Z", "AQ_total_Z")
  d$outcome <- y
  d$interaction <- "AQ×EF"
  
  slope_diff_results[[length(slope_diff_results)+1]] <- d
  
  
  # ---- AQ × DM ----
  model_4b <- lm(as.formula(paste(
    y,"~ age + sex + AQ_total_Z + EF_composit_Z + DM_composit_Z +
        AQ_total_Z:DM_composit_Z")), data=df)
  
  d <- get_slope_difference(model_4b, "DM_composit_Z", "AQ_total_Z")
  d$outcome <- y
  d$interaction <- "AQ×DM"
  
  slope_diff_results[[length(slope_diff_results)+1]] <- d
  
  
  # ---- EF × DM ----
  model_4c <- lm(as.formula(paste(
    y,"~ age + sex + AQ_total_Z + EF_composit_Z + DM_composit_Z +
        EF_composit_Z:DM_composit_Z")), data=df)
  
  d <- get_slope_difference(model_4c, "EF_composit_Z", "DM_composit_Z")
  d$outcome <- y
  d$interaction <- "EF×DM"
  
  slope_diff_results[[length(slope_diff_results)+1]] <- d
}

slope_diff_df <- do.call(rbind, slope_diff_results)

write.csv(
  slope_diff_df,
  "slope_difference_results.csv",
  row.names = FALSE
)