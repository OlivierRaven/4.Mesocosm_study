library(tidyverse)
library(readxl)
library(janitor)
library(CawthronColours)


# load data from Calum
df1 <- read_excel("./data/olliefinn data NCC.xlsx", sheet = 1) %>%
  clean_names()

# checks and fixes
head(df1)

# ok
unique(df1$treatment)

# a mess - these should be the replicate tanks
unique(df1$fish_id)

# create a "rep" variable instead of fish_id
df_list <- split(df1, interaction(df1$treatment, df1$light_regime, drop = TRUE))

df1 <- lapply(df_list, function(x){
  x$rep = rep(1:(nrow(x)/21), each = 21)
  x
}) %>%
  bind_rows() %>%
  select(-fish_id)

# ok
unique(df1$rep)

# ok
unique(df1$time_min)

# ok
unique(df1$time_zone)

# ok
unique(df1$response)

# ok
unique(df1$light_regime)


# plot the data ----------------------------------------------------------------

head(df1)

df_plot <- df1 %>%
  group_by(treatment, light_regime, time_zone, time_min, response) %>%
  summarise(tot = n()) %>%
  group_by(treatment, light_regime, time_zone, time_min) %>%
  mutate(tot2 = sum(tot),
         prop = tot/tot2) %>%
  filter(response == "REF")


rects <- data.frame(
  xmin = c(-10, 60, 240),
  xmax = c(60, 240, 320),
  ymin = -Inf,
  ymax = Inf,
  period = c("before", "during", "after"))

# plot
ggplot() +
  geom_point(data = df_plot, aes(time_min, prop, colour = treatment)) +
  geom_line(data = df_plot, aes(time_min, prop, colour = treatment)) +
  facet_wrap(~light_regime) +
  geom_rect(
    data = rects,
    aes(xmin = xmin,
        xmax = xmax,
        ymin = ymin,
        ymax = ymax,
        fill = period),
    alpha = 0.4,
    colour = NA) +
  labs(x = "Time (minutes)",
       y = "Prop. in refugia",
       colour = "Treatment",
       fill = "Period") +
  scale_fill_manual(values = RColorBrewer::brewer.pal(8, "Set2")[c(6:8)],
                    breaks = c("after", "during", "before")) +
  scale_colour_manual(values = RColorBrewer::brewer.pal(3, "Set2")) +
  scale_x_continuous(expand = c(0,0)) +
  theme_bw() +
  theme(legend.position = "bottom")

# save figure
ggsave("./outputs/plot_raw_data.png",
       width = 16,
       height = 10,
       units = "cm",
       dpi = 600)

# Fit GLMM ---------------------------------------------------------------------
# Packages
library(lme4)
library(DHARMa)
library(emmeans)
library(ggeffects)

# Ensure factors
df1$treatment    <- factor(df1$treatment)
df1$time_zone    <- factor(df1$time_zone)
df1$response     <- factor(df1$response)
df1$light_regime <- factor(df1$light_regime)
df1$rep <- factor(df1$rep)


plogis(0.5)

# set the default for each
df1$treatment <- relevel(df1$treatment, "control")
df1$time_zone <- relevel(df1$time_zone, "BEFORE")
df1$light_regime     <- relevel(df1$light_regime, "LIGHT")


# set factor ordering
df1$time_zone <- factor(df1$time_zone, levels = c("BEFORE", "DURING", "AFTER"))

# scale time so the model behaves better numerically
df1$time_sc <- scale(df1$time_min, center = TRUE, scale = TRUE)

# fit model
m_time <- glmer(
  response ~ treatment * light_regime * time_zone + time_sc +
    #  (1 + time_sc | rep), # variance for random effect time_sc was near 0.0 so dropped it.
    (1 | rep),
  data = df1,
  family = binomial,
  control = glmerControl(
    optimizer = "bobyqa",
    optCtrl = list(maxfun = 2e5)
  ))

sum_mod <- summary(m_time)
sum_mod


# extract table of fixed effects  to save --------------------------------------
fixed_df <- as.data.frame(sum_mod$coefficients)
names(fixed_df) <- c("estimate", "std_error", "z_value", "p_value")

fixed_df$term <- rownames(fixed_df)
rownames(fixed_df) <- NULL

# Reorder columns
fixed_df <- fixed_df[, c("term", "estimate", "std_error", "z_value", "p_value")]


fixed_df <- fixed_df %>%
  mutate_if(is.numeric, round, digits = 3)



write.csv(fixed_df, "./outputs/table_glmm_fixed_effects.csv", row.names = FALSE)

# 2. Basic diagnostics - all ok ------------------------------------------------
sim <- simulateResiduals(m_time)
plot(sim)
testDispersion(sim)
testZeroInflation(sim)

# 3. Post-hoc comparisons
# P value adjustment: tukey method for comparing a family of 3 estimates 

# Predator differences within each hour ----------------------------------------
emm_pred_by_hour <- emmeans(m_time, ~ treatment | time_zone * light_regime)
pairs(emm_pred_by_hour)

# convert to dataframe
emm_df <- as.data.frame(pairs(emm_pred_by_hour))

# save
write.csv(emm_df, "./outputs/table_glmm_post_hoc_predator.csv",
          row.names = F)

# Effect of hour within each predator ------------------------------------------
emm_hour_by_pred <- emmeans(m_time, ~ time_zone | treatment * light_regime)
pairs(emm_hour_by_pred)

# convert to df for saving
emm_df2 <- as.data.frame(pairs(emm_hour_by_pred))

# save
write.csv(emm_df2, "./outputs/table_glmm_post_hoc_period.csv",
          row.names = F)

# 4. Plot predicted probabilities ----------------------------------------------
pred <- ggeffect(m_time, terms = c("time_zone", "treatment", "light_regime"),
                 bias_correction = T)

# plot
ggplot(pred, aes(x = x, y = predicted,
                 colour = group, group = group)) +
  geom_line(linewidth = 1) +
  geom_point(size = 2.5) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high, fill = group),
              alpha = 0.15, colour = NA) +
  facet_wrap(~ facet, labeller = label_both) +
  scale_x_discrete("Time") +
  scale_y_continuous("Probability of hiding in refugia", limits = c(0, 1)) +
  labs(colour = "Predator", fill = "Predator") +
  theme_bw(base_size = 13) +
  theme(legend.position = "bottom")

# save figure
ggsave("./outputs/plot_glmm.png",
       width = 16,
       height = 10,
       units = "cm",
       dpi = 600)

# END --------------------------------------------------------------------------

