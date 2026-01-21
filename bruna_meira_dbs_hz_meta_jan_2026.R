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
