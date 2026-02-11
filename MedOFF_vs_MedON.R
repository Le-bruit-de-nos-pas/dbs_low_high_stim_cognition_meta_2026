
final_df <- fread("final_df_new.csv")

low_high_df <- final_df %>%
  filter(
    (condition_A == "LowHz" & condition_B == "HighHz") |
      (condition_A == "HighHz" & condition_B == "LowHz")
  )


table(low_high_df$condition_A, low_high_df$condition_B)

nrow(low_high_df)

table(low_high_df$med_state)
table(low_high_df$med_state, low_high_df$study_name)

res_med <- rma.mv(
  yi, vi,
  mods = ~ med_state,
  random = ~ 1 | study_name /  effect_id ,
  method = "REML",
  test="t",
  data = low_high_df
)

summary(res_med)


# Multivariate Meta-Analysis Model (k = 27; method: REML)
# 
#   logLik  Deviance       AIC       BIC      AICc   
#  -3.0986    6.1971   14.1971   19.0726   16.1971   
# 
# Variance Components:
# 
#             estim    sqrt  nlvls  fixed                factor 
# sigma^2.1  0.0000  0.0000      8     no            study_name 
# sigma^2.2  0.0000  0.0000     27     no  study_name/effect_id 
# 
# Test for Residual Heterogeneity:
# QE(df = 25) = 21.0827, p-val = 0.6880
# 
# Test of Moderators (coefficient 2):
# F(df1 = 1, df2 = 25) = 9.6175, p-val = 0.0047
# 
# Model Results:
# 
#              estimate      se     tval  df    pval    ci.lb   ci.ub     
# intrcpt       -0.0209  0.0734  -0.2843  25  0.7785  -0.1720  0.1302     
# med_stateOn    0.3282  0.1058   3.1012  25  0.0047   0.1102  0.5461  ** 
# 
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1


predict(res_med, newmods = model.matrix(~ med_state,
                                        data = data.frame(med_state = c("Off","On")))[, -1])

#      pred     se   ci.lb  ci.ub   pi.lb  pi.ub 
# 1 -0.0209 0.0734 -0.1720 0.1302 -0.1720 0.1302 
# 2  0.3073 0.0763  0.1503 0.4644  0.1503 0.4644 


res_rve_med <- robu(
  yi ~ med_state,
  data = low_high_df,
  studynum = effect_id,
  var.eff.size = vi,
  modelweights = "CORR",
  rho = 0.6
)

res_rve_med

# RVE: Correlated Effects Model with Small-Sample Corrections 
# 
# Model: yi ~ med_state 
# 
# Number of studies = 27 
# Number of outcomes = 27 (min = 1 , mean = 1 , median = 1 , max = 1 )
# Rho = 0.6 
# I.sq = 0 
# Tau.sq = 0 
# 
#                Estimate StdErr t-value   dfs P(|t|>) 95% CI.L 95% CI.U Sig
# 1 X.Intercept.  -0.0209 0.0837  -0.249  9.82 0.80832   -0.208    0.166    
# 2  med_stateOn   0.3282 0.1013   3.240 20.44 0.00402    0.117    0.539 ***
# ---
# Signif. codes: < .01 *** < .05 ** < .10 *
# ---
# Note: If df < 4, do not trust the results



df_off <- subset(low_high_df, med_state == "Off")
df_on  <- subset(low_high_df, med_state == "On")


res_off <- rma.mv(
  yi, vi,
  random = ~ 1 | study_name / effect_id,
  method = "REML",
  test = "t",
  data = df_off
)

res_on <- rma.mv(
  yi, vi,
  random = ~ 1 | study_name / effect_id,
  method = "REML",
  test = "t",
  data = df_on
)



xlim <- c(-1.5, 2.0)

par(mfrow = c(1, 2))


forest(
  res_off, 
  slab = paste(df_off$study_name, "-", df_off$`Cog Task`),
  xlab = "Hedges' g (Low Hz vs High Hz) \n (positive = better cognitive performance for [Low Hz])",
  main = "OFF medication",
  annotate = TRUE, addfit = TRUE,
  showweights = TRUE, header = TRUE,
  cex = 0.85,
  col = "white", border = "firebrick",
  colout = "firebrick"
)

forest(
  res_on, xlim = xlim,
  slab = paste(df_on$study_name, "-", df_on$`Cog Task`),
  xlab = "Hedges' g (Low Hz vs High Hz) \n (positive = better cognitive performance for [Low Hz])",
  main = "ON medication",
  annotate = TRUE, addfit = TRUE,
  showweights = TRUE, header = TRUE,
  cex = 0.85,
  col = "white", border = "firebrick",
  colout = "firebrick"
)

par(mfrow = c(1, 1))
# 1400 600




vlow_high_df <- final_df %>%
  filter(
    (condition_A == "VeryLowHz" & condition_B == "HighHz") |
      (condition_A == "HighHz" & condition_B == "VeryLowHz")
  )



table(vlow_high_df$condition_A, vlow_high_df$condition_B)

nrow(vlow_high_df)

table(vlow_high_df$med_state)


vlow_high_df <- vlow_high_df %>% filter(med_state!="")


table(vlow_high_df$med_state, vlow_high_df$study_name)

res_med <- rma.mv(
  yi, vi,
  mods = ~ med_state,
  random = ~ 1 | study_name /  effect_id ,
  method = "REML",
  test="t",
  data = vlow_high_df
)

summary(res_med)


# Multivariate Meta-Analysis Model (k = 57; method: REML)
# 
#   logLik  Deviance       AIC       BIC      AICc   
# -59.2784  118.5568  126.5568  134.5862  127.3568   
# 
# Variance Components:
# 
#             estim    sqrt  nlvls  fixed                factor 
# sigma^2.1  0.0000  0.0000     10     no            study_name 
# sigma^2.2  0.3499  0.5915     57     no  study_name/effect_id 
# 
# Test for Residual Heterogeneity:
# QE(df = 55) = 237.9452, p-val < .0001
# 
# Test of Moderators (coefficient 2):
# F(df1 = 1, df2 = 55) = 3.0467, p-val = 0.0865
# 
# Model Results:
# 
#              estimate      se     tval  df    pval    ci.lb   ci.ub    
# intrcpt        0.2718  0.1259   2.1584  55  0.0353   0.0194  0.5242  * 
# med_stateOn   -0.3033  0.1737  -1.7455  55  0.0865  -0.6515  0.0449  . 
# 
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

predict(res_med, newmods = model.matrix(~ med_state,
                                        data = data.frame(med_state = c("Off","On")))[, -1])

#      pred     se   ci.lb  ci.ub   pi.lb  pi.ub 
# 1  0.2718 0.1259  0.0194 0.5242 -0.9402 1.4838 
# 2 -0.0315 0.1197 -0.2713 0.2084 -1.2409 1.1780


res_rve_med <- robu(
  yi ~ med_state,
  data = vlow_high_df,
  studynum = effect_id,
  var.eff.size = vi,
  modelweights = "CORR",
  rho = 0.6
)

res_rve_med

# RVE: Correlated Effects Model with Small-Sample Corrections 
# 
# Model: yi ~ med_state 
# 
# Number of studies = 57 
# Number of outcomes = 57 (min = 1 , mean = 1 , median = 1 , max = 1 )
# Rho = 0.6 
# I.sq = 76.88544 
# Tau.sq = 0.2226238 
# 
#                Estimate StdErr t-value  dfs  P(|t|>) 95% CI.L 95% CI.U Sig
# 1 X.Intercept.    0.268 0.0656    4.08 24.8 0.000407    0.133   0.4028 ***
# 2  med_stateOn   -0.286 0.1602   -1.78 53.5 0.080270   -0.607   0.0356   *
# ---
# Signif. codes: < .01 *** < .05 ** < .10 *
# ---
# Note: If df < 4, do not trust the results



df_off <- subset(vlow_high_df, med_state == "Off")
df_on  <- subset(vlow_high_df, med_state == "On")


res_off <- rma.mv(
  yi, vi,
  random = ~ 1 | study_name / effect_id,
  method = "REML",
  test = "t",
  data = df_off
)

res_on <- rma.mv(
  yi, vi,
  random = ~ 1 | study_name / effect_id,
  method = "REML",
  test = "t",
  data = df_on
)



xlim <- c(-1.5, 2.0)

par(mfrow = c(1, 2))


forest(
  res_off, 
  slab = paste(df_off$study_name, "-", df_off$`Cog Task`),
  xlab = "Hedges' g (Very Low Hz vs High Hz) \n (positive = better cognitive performance for [Very Low Hz])",
  main = "OFF medication",
  annotate = TRUE, addfit = TRUE,
  showweights = TRUE, header = TRUE,
  cex = 0.85,
  col = "white", border = "firebrick",
  colout = "firebrick"
)

forest(
  res_on, xlim = xlim,
  slab = paste(df_on$study_name, "-", df_on$`Cog Task`),
  xlab = "Hedges' g (Very Low Hz vs High Hz) \n (positive = better cognitive performance for [Very Low Hz])",
  main = "ON medication",
  annotate = TRUE, addfit = TRUE,
  showweights = TRUE, header = TRUE,
  cex = 0.85,
  col = "white", border = "firebrick",
  colout = "firebrick"
)

par(mfrow = c(1, 1))
# 1400 600













off_high_df <- final_df %>%
  filter( (condition_A == "Off" & condition_B == "HighHz") )


table(off_high_df$condition_A, off_high_df$condition_B)

nrow(off_high_df)

table(off_high_df$med_state)


unique(off_high_df$med_state)

off_high_df <- off_high_df %>% filter(med_state!="")


table(off_high_df$med_state, off_high_df$study_name)

res_med <- rma.mv(
  yi, vi,
  mods = ~ med_state,
  random = ~ 1 | study_name /  effect_id ,
  method = "REML",
  test="t",
  data = off_high_df
)

summary(res_med)

# 
# Multivariate Meta-Analysis Model (k = 57; method: REML)
# 
#   logLik  Deviance       AIC       BIC      AICc
# -59.2784  118.5568  126.5568  134.5862  127.3568
# 
# Variance Components:
# 
#             estim    sqrt  nlvls  fixed                factor
# sigma^2.1  0.0000  0.0000     10     no            study_name
# sigma^2.2  0.3499  0.5915     57     no  study_name/effect_id
# 
# Test for Residual Heterogeneity:
# QE(df = 55) = 237.9452, p-val < .0001
# 
# Test of Moderators (coefficient 2):
# F(df1 = 1, df2 = 55) = 3.0467, p-val = 0.0865
# 
# Model Results:
# 
#              estimate      se     tval  df    pval    ci.lb   ci.ub
# intrcpt        0.2718  0.1259   2.1584  55  0.0353   0.0194  0.5242  *
# med_stateOn   -0.3033  0.1737  -1.7455  55  0.0865  -0.6515  0.0449  .
# 
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

predict(res_med, newmods = model.matrix(~ med_state,
                                        data = data.frame(med_state = c("Off","On")))[, -1])

#      pred     se   ci.lb  ci.ub   pi.lb  pi.ub 
# 1  0.0469 0.0942 -0.1417 0.2354 -0.7664 0.8601 
# 2 -0.1795 0.0955 -0.3708 0.0118 -0.9934 0.6343 


res_rve_med <- robu(
  yi ~ med_state,
  data = off_high_df,
  studynum = effect_id,
  var.eff.size = vi,
  modelweights = "CORR",
  rho = 0.6
)

res_rve_med

# RVE: Correlated Effects Model with Small-Sample Corrections 
# 
# Model: yi ~ med_state 
# 
# Number of studies = 59 
# Number of outcomes = 59 (min = 1 , mean = 1 , median = 1 , max = 1 )
# Rho = 0.6 
# I.sq = 65.4669 
# Tau.sq = 0.1208357 
# 
#                Estimate StdErr t-value  dfs P(|t|>) 95% CI.L 95% CI.U Sig
# 1 X.Intercept.   0.0528 0.0478    1.11 26.5  0.2788  -0.0453  0.15093    
# 2  med_stateOn  -0.2327 0.1209   -1.93 55.9  0.0592  -0.4749  0.00939   *
# ---
# Signif. codes: < .01 *** < .05 ** < .10 *
# ---
# Note: If df < 4, do not trust the results


df_off <- subset(off_high_df, med_state == "Off")
df_on  <- subset(off_high_df, med_state == "On")


res_off <- rma.mv(
  yi, vi,
  random = ~ 1 | study_name / effect_id,
  method = "REML",
  test = "t",
  data = df_off
)

res_on <- rma.mv(
  yi, vi,
  random = ~ 1 | study_name / effect_id,
  method = "REML",
  test = "t",
  data = df_on
)



xlim <- c(-1.5, 2.0)

par(mfrow = c(1, 2))


forest(
  res_off, 
  slab = paste(df_off$study_name, "-", df_off$`Cog Task`),
  xlab = "Hedges' g (OFF vs High Hz) \n (positive = better cognitive performance for [OFF])",
  main = "OFF medication",
  annotate = TRUE, addfit = TRUE,
  showweights = TRUE, header = TRUE,
  cex = 0.85,
  col = "white", border = "firebrick",
  colout = "firebrick"
)

forest(
  res_on, xlim = xlim,
  slab = paste(df_on$study_name, "-", df_on$`Cog Task`),
  xlab = "Hedges' g (OFF vs High Hz) \n (positive = better cognitive performance for [OFF])",
  main = "ON medication",
  annotate = TRUE, addfit = TRUE,
  showweights = TRUE, header = TRUE,
  cex = 0.85,
  col = "white", border = "firebrick",
  colout = "firebrick"
)

par(mfrow = c(1, 1))
# 1400 600





off_low_df <- final_df %>%
  filter( (condition_A == "Off" & condition_B == "LowHz") | (condition_A == "Off" & condition_B == "VeryLowHz")  )


table(off_low_df$condition_A, off_low_df$condition_B)

nrow(off_low_df)

table(off_low_df$med_state)


unique(off_low_df$med_state)

off_low_df <- off_low_df %>% filter(med_state!="")


table(off_low_df$med_state, off_low_df$study_name)

res_med <- rma.mv(
  yi, vi,
  mods = ~ med_state,
  random = ~ 1 | study_name /  effect_id ,
  method = "REML",
  test="t",
  data = off_low_df
)

summary(res_med)

# Multivariate Meta-Analysis Model (k = 61; method: REML)
# 
#   logLik  Deviance       AIC       BIC      AICc   
# -39.6513   79.3026   87.3026   95.6127   88.0433   
# 
# Variance Components:
# 
#             estim    sqrt  nlvls  fixed                factor 
# sigma^2.1  0.0121  0.1100     12     no            study_name 
# sigma^2.2  0.0717  0.2677     61     no  study_name/effect_id 
# 
# Test for Residual Heterogeneity:
# QE(df = 59) = 141.6728, p-val < .0001
# 
# Test of Moderators (coefficient 2):
# F(df1 = 1, df2 = 59) = 0.6771, p-val = 0.4139
# 
# Model Results:
# 
#              estimate      se     tval  df    pval    ci.lb    ci.ub    
# intrcpt       -0.2197  0.0830  -2.6486  59  0.0104  -0.3857  -0.0537  * 
# med_stateOn    0.1006  0.1223   0.8228  59  0.4139  -0.1441   0.3454    
# 
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

predict(res_med, newmods = model.matrix(~ med_state,
                                        data = data.frame(med_state = c("Off","On")))[, -1])

#    pred     se   ci.lb   ci.ub   pi.lb  pi.ub 
# 1 -0.2197 0.0830 -0.3857 -0.0537 -0.8222 0.3828 
# 2 -0.1191 0.0899 -0.2989  0.0608 -0.7255 0.4874 


res_rve_med <- robu(
  yi ~ med_state,
  data = off_low_df,
  studynum = effect_id,
  var.eff.size = vi,
  modelweights = "CORR",
  rho = 0.6
)

res_rve_med
# 
# RVE: Correlated Effects Model with Small-Sample Corrections 
# 
# Model: yi ~ med_state 
# 
# Number of studies = 61 
# Number of outcomes = 61 (min = 1 , mean = 1 , median = 1 , max = 1 )
# Rho = 0.6 
# I.sq = 58.35475 
# Tau.sq = 0.08881711 
# 
#                Estimate StdErr t-value  dfs  P(|t|>) 95% CI.L 95% CI.U Sig
# 1 X.Intercept.  -0.1968 0.0459   -4.29 28.3 0.000191   -0.291   -0.103 ***
# 2  med_stateOn   0.0482 0.1047    0.46 57.4 0.647201   -0.161    0.258    
# ---
# Signif. codes: < .01 *** < .05 ** < .10 *
# ---
# Note: If df < 4, do not trust the results

df_off <- subset(off_low_df, med_state == "Off")
df_on  <- subset(off_low_df, med_state == "On")


res_off <- rma.mv(
  yi, vi,
  random = ~ 1 | study_name / effect_id,
  method = "REML",
  test = "t",
  data = df_off
)

res_on <- rma.mv(
  yi, vi,
  random = ~ 1 | study_name / effect_id,
  method = "REML",
  test = "t",
  data = df_on
)



xlim <- c(-1.5, 2.0)

par(mfrow = c(1, 2))


forest(
  res_off, 
  slab = paste(df_off$study_name, "-", df_off$`Cog Task`),
  xlab = "Hedges' g (OFF vs Low Hz) \n (positive = better cognitive performance for [OFF])",
  main = "OFF medication",
  annotate = TRUE, addfit = TRUE,
  showweights = TRUE, header = TRUE,
  cex = 0.85,
  col = "white", border = "firebrick",
  colout = "firebrick"
)

forest(
  res_on, xlim = xlim,
  slab = paste(df_on$study_name, "-", df_on$`Cog Task`),
  xlab = "Hedges' g (OFF vs Low Hz) \n (positive = better cognitive performance for [OFF])",
  main = "ON medication",
  annotate = TRUE, addfit = TRUE,
  showweights = TRUE, header = TRUE,
  cex = 0.85,
  col = "white", border = "firebrick",
  colout = "firebrick"
)

par(mfrow = c(1, 1))
# 1400 600



