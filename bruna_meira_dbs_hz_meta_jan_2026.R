library(data.table)
library(tidyverse)
library(readxl)
library(purrr)
library(metafor)


# Data pre-processing -------------------------

final_df <- read_excel("LFS.HFS.toanalyse.new.xlsx", sheet = "Folha1")

names(final_df)

length(unique(final_df$study_name)) # 22

final_df %>% select(study_name, design) %>% distinct() %>% count() # 23

unique(final_df$contrast)

assumed_r <- 0.5 

unique(final_df$design)

final_df <- final_df %>%
  mutate(
    design_type = if_else(
      grepl("Within", design, ignore.case = TRUE),
      "within",
      "between"
    )
  )

# Within-subject Hedges g (small sample correct SMD)

hedges_g_within <- function(m1, m2, sd1, sd2, n, r) {
  sd_diff <- sqrt(sd1^2 + sd2^2 - 2*r*sd1*sd2)
  dz <- (m1 - m2) / sd_diff
  J <- 1 - (3 / (4*n - 1))
  g <- dz * J
  var_g <- (1/n) + (g^2 / (2*n))
  list(g = g, v = var_g)
}


# Between-subject Hedges g
hedges_g_between <- function(m1, m2, sd1, sd2, n1, n2) {
  s_pooled <- sqrt(
    ((n1 - 1)*sd1^2 + (n2 - 1)*sd2^2) / (n1 + n2 - 2)
  )
  d <- (m1 - m2) / s_pooled
  J <- 1 - (3 / (4*(n1 + n2) - 9))
  g <- d * J
  var_g <- (n1 + n2)/(n1*n2) + (g^2/(2*(n1 + n2)))
  list(g = g, v = var_g)
}



final_df <- final_df %>%
  mutate(
    mean_A = as.numeric(mean_A),
    mean_B = as.numeric(mean_B),
    sd_A   = as.numeric(sd_A),
    sd_B   = as.numeric(sd_B),
    n_A    = as.numeric(n_A),
    n_B    = as.numeric(n_B)
  )


final_df <- final_df %>%
  rowwise() %>%
  mutate(
    es = list(
      if (design_type == "within") {
        hedges_g_within(
          mean_A, mean_B,
          sd_A, sd_B,
          n_A,
          assumed_r
        )
      } else {
        hedges_g_between(
          mean_A, mean_B,
          sd_A, sd_B,
          n_A, n_B
        )
      }
    ),
    yi = es$g,
    vi = es$v
  ) %>%
  ungroup() %>%
  select(-es)



# SMD/hedges_g sign correction


final_df <- final_df %>%
  mutate(
    yi = if_else(higher_is_better=="TRUE", yi, -yi)
  )


summary(final_df$yi)
summary(final_df$vi)


# How to treat the nested structure
# Use study_id to treat all contrasts as independendt, study_name as dependent

final_df <- final_df %>%
  mutate(
    study_id = as.factor(study_id),
    effect_id = row_number()
  )


# final_df <- final_df %>%
#   group_by(study_name) %>%         # group by study
#   mutate(effect_id = row_number()) %>%  # effect_id unique within study
#   ungroup() %>%
#   mutate(study_name = as.factor(study_name),
#          effect_id = as.factor(effect_id))


fwrite(final_df, "final_df_new.csv")


# ----------------

# LOW vs HIGH Overall and Main Domains  ------------


low_high_df <- final_df %>%
  filter(
    (condition_A == "LowHz" & condition_B == "HighHz") |
      (condition_A == "HighHz" & condition_B == "LowHz")
  )


table(low_high_df$condition_A, low_high_df$condition_B)

nrow(low_high_df)


res_lh <- rma.mv(
  yi, vi,
  random = ~ 1 | study_name / study_id ,
  method = "REML",
  test="t",
  data = low_high_df
)

summary(res_lh)


tau2 <- res_lh$sigma2
I2 <- 100 * tau2 / (tau2 + mean(low_high_df$vi))

tau2
I2


run_with_r_lh <- function(r_val) {
  low_high_df %>%
    rowwise() %>%
    mutate(
      es = list(
        if (design_type == "within") {
          hedges_g_within(mean_A, mean_B, sd_A, sd_B, n_A, r_val)
        } else {
          hedges_g_between(mean_A, mean_B, sd_A, sd_B, n_A, n_B)
        }
      ),
      yi = es$g,
      vi = es$v,
      yi = if_else(higher_is_better=="TRUE", yi, -yi)  # flip sign
    ) %>%
    ungroup() %>%
    rma.mv(
      yi, vi,
      random = ~ 1 |study_name,
      method = "REML",
      data = .
    )
}

summary(run_with_r_lh(0.3)) #012
summary(run_with_r_lh(0.5)) #0.14
summary(run_with_r_lh(0.7)) # 0.17


meta::forest(
  res_lh,
  annotate=TRUE, addfit=TRUE, 
  showweights=TRUE, header=TRUE,
  slab = paste(low_high_df$study_name, "-", low_high_df$`Cog Task`),
  xlab = "Hedges' g (positive = better cognitive performance for [Low Hz])",
  cex = 0.9,
  col="white", border="firebrick",
  colout="firebrick"
)




# Define a sequence of r values to test
r_vals <- seq(0, 0.9, by = 0.1)


# Function to compute meta-analysis for a given r
run_with_r <- function(r_val) {
  low_high_df %>%
    rowwise() %>%
    mutate(
      es = list(
        if (design_type == "within") {
          hedges_g_within(mean_A, mean_B, sd_A, sd_B, n_A, r_val)
        } else {
          hedges_g_between(mean_A, mean_B, sd_A, sd_B, n_A, n_B)
        }
      ),
      yi = es$g,
      vi = es$v,
      yi = if_else(higher_is_better=="TRUE", yi, -yi)   # flip sign if needed
    ) %>%
    ungroup() %>%
    rma.mv(
      yi, vi,
      random = ~ 1 | study_id / effect_id,
      method = "REML",
      data = .
    ) %>%
    { tibble(
      r = r_val,
      estimate = .$beta[1],
      se = .$se,
      ci.lb = .$ci.lb,
      ci.ub = .$ci.ub
    ) }
}

# Run the sensitivity analysis across all r values
results <- map_df(r_vals, run_with_r)

# Plot effect size vs r
results %>%
  ggplot(aes(x = r, y = estimate)) +
  ylim(-0.5,0.5) +
  geom_line(color = "deepskyblue4", size = 1) +
  geom_ribbon(aes(ymin = ci.lb, ymax = ci.ub), alpha = 0.1, fill = "deepskyblue4") +
  geom_hline(yintercept = 0, linetype = "dashed", colour="firebrick", size=1) +
  labs(
    x = "\n Assumed correlation r",
    y = "Effect size (Hedges' g) \n",
    title = "Low Hz vs High Hz \nSensitivity of effect size to assumed correlation r"
  ) +
  theme(axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        legend.position = "none") +
  theme(panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        strip.text = element_blank(),
        axis.line = element_blank(),
        axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 10),
        axis.title.x = element_text(size = 10, vjust = -0.5),
        axis.title.y = element_text(size = 10, vjust = -0.5),
        plot.margin = margin(5, 5, 5, 5, "pt")) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))




table(low_high_df$Main_cog_domain)


# Identify domains with at least 2 effect sizes
domain_counts <- table(low_high_df$Main_cog_domain)
valid_domains <- names(domain_counts[domain_counts >= 2])

# Create empty list to store results
domain_results <- list()
forest_plots <- list()

for(d in valid_domains){
  
  df <- low_high_df %>% filter(Main_cog_domain == d)
  
  # Run multilevel meta-analysis
  res <- rma.mv(
    yi = yi,
    V = vi,
    random = ~1 | study_name,
    method = "REML",
    data = df
  )
  
  # Store results
  domain_results[[d]] <- tibble(
    domain = d,
    estimate = res$beta[1],
    se = res$se,
    ci.lb = res$ci.lb,
    ci.ub = res$ci.ub,
    pval = res$pval
  )
  
  # Create forest plot
  forest(res,
         slab = paste(df$study_name, "-", df$`Cog Task`),
         xlab = "Hedges' g (positive = better cognitive performance for [Low Hz])",
         main = paste("\n Forest plot:", d),
         annotate=TRUE, addfit=TRUE, 
         showweights=TRUE, header=TRUE,
         cex = 0.9,
         col="white", border="firebrick",
         colout="firebrick")
  
  forest_plots[[d]] <- res
}

# Combine all results
domain_results_df <- bind_rows(domain_results)
domain_results_df


library(robumeta)


# Small-sample corrected RVE
res_rve <- robu(
  formula = yi ~ 1,           # intercept-only model
  data = low_high_df,
  studynum = study_id,        # clusters by study
  var.eff.size = vi,          # variance of effect sizes
  modelweights = "CORR",      # correlated effects model
  rho = 0.6                   # assumed correlation between effect sizes within study
)

summary(res_rve)
res_rve




# LEAVE ON OUT CV 

studies <- unique(low_high_df$study_name)

loso_results <- lapply(studies, function(s) {
  
  res <- rma.mv(
    yi, vi,
    random = ~ 1 | study_id / effect_id,
    method = "REML",
    data = subset(low_high_df, study_name != s)
  )
  
  tibble(
    left_out = s,
    estimate = res$beta[1],
    se = res$se,
    ci.lb = res$ci.lb,
    ci.ub = res$ci.ub
  )
})

loso_df <- bind_rows(loso_results)
loso_df

library(metafor)

forest(
  x = loso_df$estimate,
  sei = loso_df$se,
  slab = paste("Leave out:", loso_df$left_out),
  xlab = "Pooled Hedges' g (Leave-One-Out)",
  alim = c(-0.2, 0.4),
  at = seq(-0.2, 0.4, 0.1),
  cex = 0.8
)




# -----------------------------
# VERY LOW VS VERY HIGH Overall and Main Domains-----------



vlow_high_df <- final_df %>%
  filter(
    (condition_A == "VeryLowHz" & condition_B == "HighHz") |
      (condition_A == "HighHz" & condition_B == "VeryLowHz")
  )


table(vlow_high_df$condition_A, vlow_high_df$condition_B)

nrow(vlow_high_df)

unique(vlow_high_df$study_name)


res_vlh <- rma.mv(
  yi, vi,
  random = ~ 1 | study_id / effect_id,
  method = "REML",
  data = vlow_high_df
)

summary(res_vlh)




tau2 <- res_vlh$sigma2
I2 <- 100 * tau2 / (tau2 + mean(vlow_high_df$vi))
tau2
I2



run_with_r_vlh <- function(r_val) {
  vlow_high_df %>%
    rowwise() %>%
    mutate(
      es = list(
        if (design_type == "within") {
          hedges_g_within(mean_A, mean_B, sd_A, sd_B, n_A, r_val)
        } else {
          hedges_g_between(mean_A, mean_B, sd_A, sd_B, n_A, n_B)
        }
      ),
      yi = es$g,
      vi = es$v,
      yi = if_else(higher_is_better=="TRUE", yi, -yi)  # flip sign
    ) %>%
    ungroup() %>%
    rma.mv(
      yi, vi,
      random = ~ 1 | study_name,
      method = "REML",
      data = .
    )
}

summary(run_with_r_vlh(0.3))
summary(run_with_r_vlh(0.5))

summary(run_with_r_vlh(0.7))



meta::forest(
  res_vlh,
  annotate=TRUE, addfit=TRUE, 
  showweights=TRUE, header=TRUE,
  slab = paste(vlow_high_df$study_name, "-", vlow_high_df$`Cog Task`),
  xlab = "Hedges' g (positive = better cognitive performance for [Very Low Hz])",
  cex = 0.9,
  col="white", border="firebrick",
  colout="firebrick"
)



# Define a sequence of r values to test
r_vals <- seq(0, 0.9, by = 0.1)


# Function to compute meta-analysis for a given r
run_with_r <- function(r_val) {
  vlow_high_df %>%
    rowwise() %>%
    mutate(
      es = list(
        if (design_type == "within") {
          hedges_g_within(mean_A, mean_B, sd_A, sd_B, n_A, r_val)
        } else {
          hedges_g_between(mean_A, mean_B, sd_A, sd_B, n_A, n_B)
        }
      ),
      yi = es$g,
      vi = es$v,
      yi = if_else(higher_is_better=="TRUE", yi, -yi)   # flip sign if needed
    ) %>%
    ungroup() %>%
    rma.mv(
      yi, vi,
      random = ~ 1 | study_name,
      method = "REML",
      data = .
    ) %>%
    { tibble(
      r = r_val,
      estimate = .$beta[1],
      se = .$se,
      ci.lb = .$ci.lb,
      ci.ub = .$ci.ub
    ) }
}

# Run the sensitivity analysis across all r values
results <- map_df(r_vals, run_with_r)

# Plot effect size vs r
results %>%
  ggplot(aes(x = r, y = estimate)) +
  ylim(-0.5,1) +
  geom_line(color = "deepskyblue4", size = 1) +
  geom_ribbon(aes(ymin = ci.lb, ymax = ci.ub), alpha = 0.1, fill = "deepskyblue4") +
  geom_hline(yintercept = 0, linetype = "dashed", colour="black", size=1) +
  labs(
    x = "\n Assumed correlation r",
    y = "Effect size (Hedges' g) \n",
    title = "Very Low Hz vs High Hz \nSensitivity of effect size to assumed correlation r"
  ) +
  theme(axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        legend.position = "none") +
  theme(panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        strip.text = element_blank(),
        axis.line = element_blank(),
        axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 10),
        axis.title.x = element_text(size = 10, vjust = -0.5),
        axis.title.y = element_text(size = 10, vjust = -0.5),
        plot.margin = margin(5, 5, 5, 5, "pt")) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


table(vlow_high_df$Main_cog_domain)


# Identify domains with at least 2 effect sizes
domain_counts <- table(vlow_high_df$Main_cog_domain)
valid_domains <- names(domain_counts[domain_counts >= 2])

# Create empty list to store results
domain_results <- list()
forest_plots <- list()

for(d in valid_domains){
  
  df <- vlow_high_df %>% filter(Main_cog_domain == d)
  
  # Run multilevel meta-analysis
  res <- rma.mv(
    yi = yi,
    V = vi,
    random = ~1 | study_name,
    method = "REML",
    data = df
  )
  
  # Store results
  domain_results[[d]] <- tibble(
    domain = d,
    estimate = res$beta[1],
    se = res$se,
    ci.lb = res$ci.lb,
    ci.ub = res$ci.ub,
    pval = res$pval
  )
  
  # Create forest plot
  forest(res,
         slab = paste(df$study_name, "-", df$`Cog Task`),
         xlab = "Hedges' g (positive = better cognitive performance for [Very Low Hz])",
         main = paste("\n Forest plot:", d),
         annotate=TRUE, addfit=TRUE, 
         showweights=TRUE, header=TRUE,
         cex = 0.9,
         col="white", border="firebrick",
         colout="firebrick")
  
  
  forest_plots[[d]] <- res
}

# Combine all results
domain_results_df <- bind_rows(domain_results)
domain_results_df



library(robumeta)


# Small-sample corrected RVE
res_rve <- robu(
  formula = yi ~ 1,           # intercept-only model
  data = vlow_high_df,
  studynum = study_id,        # clusters by study
  var.eff.size = vi,          # variance of effect sizes
  modelweights = "CORR",      # correlated effects model
  rho = 0.6                   # assumed correlation between effect sizes within study
)

summary(res_rve)
res_rve




# LEAVE ON OUT CV 

studies <- unique(vlow_high_df$study_name)

loso_results <- lapply(studies, function(s) {
  
  res <- rma.mv(
    yi, vi,
    random = ~ 1 | study_id / effect_id,
    method = "REML",
    data = subset(vlow_high_df, study_name != s)
  )
  
  tibble(
    left_out = s,
    estimate = res$beta[1],
    se = res$se,
    ci.lb = res$ci.lb,
    ci.ub = res$ci.ub
  )
})

loso_df <- bind_rows(loso_results)
loso_df




library(metafor)

forest(
  x = loso_df$estimate,
  sei = loso_df$se,
  slab = paste("Leave out:", loso_df$left_out),
  xlab = "Pooled Hedges' g (Leave-One-Out)",
  alim = c(-0.2, 0.4),
  at = seq(-0.2, 0.4, 0.1),
  cex = 0.8
)


# --------

# Other (RVE, etc.) ------------

library(metafor)
library(clubSandwich)

# Fit models
naive <- rma(yi, vi, data = vlow_high_df)
mlma <- rma.mv(yi, vi, random = ~ 1 | study_name, data = vlow_high_df)
rve_result <- robust(naive, cluster = vlow_high_df$study_name)


#df = low_high_df[vlow_high_df$cog_domain=="Verbal Fluency",]

naive <- rma.mv(
  yi = yi,
  V = vi,
  method = "REML",
  data = vlow_high_df
)

mlma <- rma.mv(yi, vi, random = ~ 1 | study_name, data = vlow_high_df)
rve_result <- robust(naive, cluster = vlow_high_df$study_name)




# Compare
comparison <- data.frame(
  Method = c("Naive", "MLMA", "RVE"),
  Coefficient = c(coef(naive), coef(mlma), coef(rve_result)),
  SE = c(naive$se, sqrt(diag(vcov(mlma)))[1], rve_result$se),
  CI_lower = c(naive$ci.lb, mlma$ci.lb, rve_result$ci.lb),
  CI_upper = c(naive$ci.ub, mlma$ci.ub, rve_result$ci.ub),
  p_value = c(naive$pval, mlma$pval, rve_result$pval)
)

print(comparison)

# ---------
