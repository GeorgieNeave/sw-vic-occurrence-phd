# Run stacked single season occupancy models in ubms
# plot cumulative detection probabilities, EVC group and binary effects of baiting 
# Matt Rees

# data formatting
library(dplyr)
library(forcats)
library(camtrapR)
library(unmarked)
# models
library(ubms)
library(rstan)
# plots
library(ggplot2)
library(ggpubr)
library(patchwork)


# Set the default theme for ggplot objects 
theme_set(theme_bw())
theme_update(panel.grid = element_blank())

# set number of cores to parallelise chains
options(mc.cores = parallel::detectCores())


# PREPARE DATA ------------------------------------------------------------
## Load records 
records <- read.csv("raw_data/all_records_merged_veg_fire_raster_dist_foxcontrol_rainfall.csv")
# drop zoi's records - too much spatial autocorrelation
records <- filter(records, data_source != "zoi")
# exclude stations left out for less than 10 days
records <- filter(records, survey_duration >= 10)

## Adjust variables
# as we only had three "Riverine Grassy Woodlands or Forests" unique survey sites --> drop em
records <- filter(records, XGROUPNAME != "Riverine Grassy Woodlands or Forests")
# as we only had 20 "Rainforests" unique survey sites, all of which interspersed in "Wet Forests" in the Otways, reclassify as "Wet Forests". 
records$XGROUPNAME <- if_else(records$XGROUPNAME == "Rainforests", "Wet or Damp Forests", as.character(records$XGROUPNAME))
# abbreviate EVC group names for plotting
records$XGROUPNAME <- if_else(records$XGROUPNAME == "Wet or Damp Forests", "Wet Forests", records$XGROUPNAME)
records$XGROUPNAME <- if_else(records$XGROUPNAME == "Riparian Scrubs or Swampy Scrubs and Woodlands", "Swampy Scrubs", records$XGROUPNAME)
records$XGROUPNAME <- substr(records$XGROUPNAME, 1, nchar(records$XGROUPNAME) - 1)
table(records$XGROUPNAME)
# make a binary baiting variable
records$bait01 <- if_else(records$foxbaits == 0, "unbaited", "baited")
# make an interaction term for region / bait01
records$region_bait01 <- paste0(records$region, "_", records$bait01)

# make camdata file (records contain all detections of species / people-set-ups / falsetriggers)
camdata <- distinct(records, station_year, .keep_all = TRUE)


# drop and transform variables
camdata <- transmute(camdata,
                     station = as.factor(station),
                     station_year = as.factor(station_year),
                     survey_duration,
                     date_start,
                     date_end,
                     region = as.factor(region),
                     bait01 = as.factor(bait01),
                     XGROUPNAME = as.factor(XGROUPNAME),
                     region_bait01 = as.factor(region_bait01))

# build a camera operation matrix using camtrapR
camop <- cameraOperation(CTtable = camdata, stationCol = "station_year", setupCol = "date_start", retrievalCol = "date_end", writecsv = FALSE, hasProblems = FALSE, dateFormat = "%Y-%m-%d")

# build unmarked detection histories (using camtrapR)
dethist_fox <- as.data.frame(detectionHistory(recordTable = records, camOp = camop, stationCol = "station_year", speciesCol = "species", recordDateTimeCol = "date_time", species = "fox", occasionLength = 1, minActiveDaysPerOccasion = 1, occasionStartTime = 0, day1 = "station", datesAsOccasionNames = FALSE, includeEffort = FALSE, scaleEffort = FALSE, timeZone = "Australia/Brisbane", writecsv = FALSE))
dethist_cat <- as.data.frame(detectionHistory(recordTable = records, camOp = camop, stationCol = "station_year", speciesCol = "species", recordDateTimeCol = "date_time", species = "cat", occasionLength = 1, minActiveDaysPerOccasion = 1, occasionStartTime = 0, day1 = "station", datesAsOccasionNames = FALSE, includeEffort = FALSE, scaleEffort = FALSE, timeZone = "Australia/Brisbane", writecsv = FALSE))
dethist_sbb <- as.data.frame(detectionHistory(recordTable = records, camOp = camop, stationCol = "station_year", speciesCol = "species", recordDateTimeCol = "date_time", species = "bandicoot_sb", occasionLength = 1, minActiveDaysPerOccasion = 1, occasionStartTime = 0, day1 = "station", datesAsOccasionNames = FALSE, includeEffort = FALSE, scaleEffort = FALSE, timeZone = "Australia/Brisbane", writecsv = FALSE))
dethist_lnp <- as.data.frame(detectionHistory(recordTable = records, camOp = camop, stationCol = "station_year", speciesCol = "species", recordDateTimeCol = "date_time", species = "potoroo_ln", occasionLength = 1, minActiveDaysPerOccasion = 1, occasionStartTime = 0, day1 = "station", datesAsOccasionNames = FALSE, includeEffort = FALSE, scaleEffort = FALSE, timeZone = "Australia/Brisbane", writecsv = FALSE))

# build unmarked frames
umf_fox = unmarkedFrameOccu(dethist_fox, siteCovs = camdata, obsCovs=NULL)
umf_cat = unmarkedFrameOccu(dethist_cat, siteCovs = camdata, obsCovs=NULL)
umf_sbb = unmarkedFrameOccu(dethist_sbb, siteCovs = camdata, obsCovs=NULL)
umf_lnp = unmarkedFrameOccu(dethist_lnp, siteCovs = camdata, obsCovs=NULL)


# FIT MODELS --------------------------------------------------------------

# fit single season occu models (random effect for repeat sampling)
fit_occu_fox <- stan_occu(~region_bait01  ~region_bait01 + XGROUPNAME + (1|station), data = umf_fox, chains = 4, iter = 10000)
fit_occu_cat <- stan_occu(~region_bait01  ~region_bait01 + XGROUPNAME + (1|station), data = umf_cat, chains = 4, iter = 10000)
fit_occu_sbb <- stan_occu(~region_bait01  ~region_bait01 + XGROUPNAME + (1|station), data = umf_sbb, chains = 4, iter = 10000)
fit_occu_lnp <- stan_occu(~region_bait01  ~region_bait01 + XGROUPNAME + (1|station), data = umf_lnp, chains = 4, iter = 10000)

# save
saveRDS(fit_occu_fox, "derived_data/fit_occu_fox.RData")
saveRDS(fit_occu_cat, "derived_data/fit_occu_cat.RData")
saveRDS(fit_occu_sbb, "derived_data/fit_occu_sbb.RData")
saveRDS(fit_occu_lnp, "derived_data/fit_occu_lnp.RData")

fit_occu_fox <- readRDS("/Users/mrees2/Dropbox/personal/matt/github/unlinked/occurrence-gams-fox-cat-sbb-lnp/derived_data/fit_occu_fox.RData")
fit_occu_cat <- readRDS("/Users/mrees2/Dropbox/personal/matt/github/unlinked/occurrence-gams-fox-cat-sbb-lnp/derived_data/fit_occu_cat.RData")
fit_occu_sbb <- readRDS("/Users/mrees2/Dropbox/personal/matt/github/unlinked/occurrence-gams-fox-cat-sbb-lnp/derived_data/fit_occu_sbb.RData")
fit_occu_lnp <- readRDS("/Users/mrees2/Dropbox/personal/matt/github/unlinked/occurrence-gams-fox-cat-sbb-lnp/derived_data/fit_occu_lnp.RData")

# plot veg effects
plot_marginal(fit_occu_fox, "state") 
plot_marginal(fit_occu_fox, "det") 
# + summaries, models checks etc. 


# PLOT: CUMULATIVE DETECTABILITY -----------------------------------------------------------

# predict detection predictions (1 day occasion)
fox_det_g_b  <- predict(fit_occu_fox, submodel = "det", newdata = expand.grid(region_bait01 = "glenelg_baited",   XGROUPNAME = "Heathy Woodland"), re.form = NA)
cat_det_g_b  <- predict(fit_occu_cat, submodel = "det", newdata = expand.grid(region_bait01 = "glenelg_baited",   XGROUPNAME = "Heathy Woodland"), re.form = NA)
sbb_det_g_b  <- predict(fit_occu_sbb, submodel = "det", newdata = expand.grid(region_bait01 = "glenelg_baited",   XGROUPNAME = "Heathy Woodland"), re.form = NA)
lnp_det_g_b  <- predict(fit_occu_lnp, submodel = "det", newdata = expand.grid(region_bait01 = "glenelg_baited",   XGROUPNAME = "Heathy Woodland"), re.form = NA)
fox_det_g_ub <- predict(fit_occu_fox, submodel = "det", newdata = expand.grid(region_bait01 = "glenelg_unbaited", XGROUPNAME = "Heathy Woodland"), re.form = NA)
cat_det_g_ub <- predict(fit_occu_cat, submodel = "det", newdata = expand.grid(region_bait01 = "glenelg_unbaited", XGROUPNAME = "Heathy Woodland"), re.form = NA)
sbb_det_g_ub <- predict(fit_occu_sbb, submodel = "det", newdata = expand.grid(region_bait01 = "glenelg_unbaited", XGROUPNAME = "Heathy Woodland"), re.form = NA)
lnp_det_g_ub <- predict(fit_occu_lnp, submodel = "det", newdata = expand.grid(region_bait01 = "glenelg_unbaited", XGROUPNAME = "Heathy Woodland"), re.form = NA)
fox_det_o_b  <- predict(fit_occu_fox, submodel = "det", newdata = expand.grid(region_bait01 = "otways_baited",    XGROUPNAME = "Heathy Woodland"), re.form = NA)
cat_det_o_b  <- predict(fit_occu_cat, submodel = "det", newdata = expand.grid(region_bait01 = "otways_baited",    XGROUPNAME = "Heathy Woodland"), re.form = NA)
sbb_det_o_b  <- predict(fit_occu_sbb, submodel = "det", newdata = expand.grid(region_bait01 = "otways_baited",    XGROUPNAME = "Heathy Woodland"), re.form = NA)
lnp_det_o_b  <- predict(fit_occu_lnp, submodel = "det", newdata = expand.grid(region_bait01 = "otways_baited",    XGROUPNAME = "Heathy Woodland"), re.form = NA)
fox_det_o_ub <- predict(fit_occu_fox, submodel = "det", newdata = expand.grid(region_bait01 = "otways_unbaited",  XGROUPNAME = "Heathy Woodland"), re.form = NA)
cat_det_o_ub <- predict(fit_occu_cat, submodel = "det", newdata = expand.grid(region_bait01 = "otways_unbaited",  XGROUPNAME = "Heathy Woodland"), re.form = NA)
sbb_det_o_ub <- predict(fit_occu_sbb, submodel = "det", newdata = expand.grid(region_bait01 = "otways_unbaited",  XGROUPNAME = "Heathy Woodland"), re.form = NA)
lnp_det_o_ub <- predict(fit_occu_lnp, submodel = "det", newdata = expand.grid(region_bait01 = "otways_unbaited",  XGROUPNAME = "Heathy Woodland"), re.form = NA)

## Estimate cumulative detection probabilities w/ 95% CI's - add to single dataframe
# function to make a df per species
cum_det <- function(det_est, species, region, treatment){
  n = 1:max(records$survey_duration)
  det_prob <- 1-(1-det_est$Predicted)^n
  det_prob_lcl <- 1-(1-det_est$'2.5%')^n 
  det_prob_ucl <- 1-(1-det_est$'97.5%')^n 
  df <- as.data.frame(cbind(n, det_prob, det_prob_lcl, det_prob_ucl))
  df$Species <- species
  df$Region <- region
  df$Region <- region
  df$Treatment <- treatment
  return(df)
}
# use the function
fox_det_g_b_cdp  <- cum_det(fox_det_g_b , "Red fox",                  "Glenelg", "Baited")
cat_det_g_b_cdp  <- cum_det(cat_det_g_b , "Feral cat",                "Glenelg", "Baited")
sbb_det_g_b_cdp  <- cum_det(sbb_det_g_b , "Southern brown bandicoot", "Glenelg", "Baited")
lnp_det_g_b_cdp  <- cum_det(lnp_det_g_b , "Long-nosed potoroo",       "Glenelg", "Baited")
fox_det_g_ub_cdp <- cum_det(fox_det_g_ub, "Red fox",                  "Glenelg", "Unbaited")
cat_det_g_ub_cdp <- cum_det(cat_det_g_ub, "Feral cat",                "Glenelg", "Unbaited")
sbb_det_g_ub_cdp <- cum_det(sbb_det_g_ub, "Southern brown bandicoot", "Glenelg", "Unbaited")
lnp_det_g_ub_cdp <- cum_det(lnp_det_g_ub, "Long-nosed potoroo",       "Glenelg", "Unbaited")
fox_det_o_b_cdp  <- cum_det(fox_det_o_b , "Red fox",                  "Otways", "Baited")
cat_det_o_b_cdp  <- cum_det(cat_det_o_b , "Feral cat",                "Otways", "Baited")
sbb_det_o_b_cdp  <- cum_det(sbb_det_o_b , "Southern brown bandicoot", "Otways", "Baited")
lnp_det_o_b_cdp  <- cum_det(lnp_det_o_b , "Long-nosed potoroo",       "Otways", "Baited")
fox_det_o_ub_cdp <- cum_det(fox_det_o_ub, "Red fox",                  "Otways", "Unbaited")
cat_det_o_ub_cdp <- cum_det(cat_det_o_ub, "Feral cat",                "Otways", "Unbaited")
sbb_det_o_ub_cdp <- cum_det(sbb_det_o_ub, "Southern brown bandicoot", "Otways", "Unbaited")
lnp_det_o_ub_cdp <- cum_det(lnp_det_o_ub, "Long-nosed potoroo",       "Otways", "Unbaited")

# combine dataframes
cd_probs <- bind_rows(fox_det_g_b_cdp, 
                      cat_det_g_b_cdp ,
                      sbb_det_g_b_cdp ,
                      lnp_det_g_b_cdp ,
                      fox_det_g_ub_cdp,
                      cat_det_g_ub_cdp,
                      sbb_det_g_ub_cdp,
                      lnp_det_g_ub_cdp,
                      fox_det_o_b_cdp ,
                      cat_det_o_b_cdp ,
                      sbb_det_o_b_cdp ,
                      lnp_det_o_b_cdp ,
                      fox_det_o_ub_cdp,
                      cat_det_o_ub_cdp,
                      sbb_det_o_ub_cdp,
                      lnp_det_o_ub_cdp)


# summaries
# average:
#cd_probs_mean <- filter(cd_probs, n == 47) %>% 
#  arrange(det_prob)
## days for 95% prob
#cd_probs_95 <- filter(cd_probs, det_prob >= 0.95 & Treatment == "Unbaited") %>% 
#  arrange(n) %>% 
#  group_by(Species, Region) %>%
#  slice_head()
#cd_probs_95

# reorder species
cd_probs$Species <- fct_relevel(cd_probs$Species, "Red fox", "Feral cat", "Southern brown bandicoot", "Long-nosed potoroo")

# plot
plot_g <- ggplot(data=filter(cd_probs, Region == "Glenelg"), aes(x=n, y=det_prob, colour = Treatment)) +
  geom_ribbon(aes(ymin=det_prob_lcl, ymax=det_prob_ucl, fill = Treatment), alpha=0.2, linetype=0) +
  ylim(0,1) +
  scale_fill_manual(values = c("tomato", "dodgerblue")) +
  scale_color_manual(values = c("tomato", "dodgerblue")) +   
  facet_wrap(~Species) +
  labs(title = "Glenelg region",x = "Survey duration (days)", y = "Pr(detection)") +
  geom_line(aes(y=det_prob, x=n), lwd = 0.8) +
  geom_vline(xintercept = mean(camdata$survey_duration), colour = "black", size = 0.6) + 
  geom_vline(xintercept = quantile(camdata$survey_duration, probs = c(0.25, 0.75)), colour = "black", size = 0.6, linetype="dotted") + 
  theme(plot.title = element_text(size=11),
        axis.title = element_text(size = 10))

plot_o <- ggplot(data=filter(cd_probs, Region == "Otways"), aes(x=n, y=det_prob, colour = Treatment)) +
  geom_ribbon(aes(ymin=det_prob_lcl, ymax=det_prob_ucl, fill = Treatment), alpha=0.2, linetype=0) +
  ylim(0,1) +
  scale_fill_manual(values = c("tomato", "dodgerblue")) +
  scale_color_manual(values = c("tomato", "dodgerblue")) +   
  facet_wrap(~Species) +
  labs(title = "Otway Ranges", x = "Survey duration (days)", y = "Pr(detection)") +
  geom_line(aes(y=det_prob, x=n), lwd = 0.8) +
  geom_vline(xintercept = mean(camdata$survey_duration), colour = "black", size = 0.6) + 
  geom_vline(xintercept = quantile(camdata$survey_duration, probs = c(0.25, 0.75)), colour = "black", size = 0.6, linetype="dotted") + 
  theme(plot.title = element_text(size=11),
        axis.title = element_text(size = 10))

plot_det <- ggarrange(plot_g, plot_o, ncol=1, nrow=2, common.legend = TRUE, legend="bottom", labels = c("a", "b"),   font.label = list(size = 13, color = "black", face = "plain", family = NULL))
plot_det

# plot / save
png("figs/detectability.png", width = 8, height = 10, res = 600, units = "in")
plot_det
dev.off()


# PLOT: DETECTABILITY -----------------------------------------------------
# make a dataframe to predict into
df <- expand.grid(region_bait01 = unique(records$region_bait01), XGROUPNAME = "Heathy Woodland")
# predict into
fox_det_pred <- predict(fit_occu_fox, newdata = df, submodel = "det", re.form = NA) %>% bind_cols(df)
cat_det_pred <- predict(fit_occu_cat, newdata = df, submodel = "det", re.form = NA) %>% bind_cols(df)
sbb_det_pred <- predict(fit_occu_sbb, newdata = df, submodel = "det", re.form = NA) %>% bind_cols(df)
lnp_det_pred <- predict(fit_occu_lnp, newdata = df, submodel = "det", re.form = NA) %>% bind_cols(df)

# combine df's
det_preds <- bind_rows(fox_det_pred, cat_det_pred, sbb_det_pred, lnp_det_pred) %>%
  mutate(Species = c(rep(("Red fox"),4), rep(("Feral cat"),4), rep(("Southern brown bandicoot"),4),  rep(("Long-nosed potoroo"),4)),
         Region = rep(c(rep(("Glenelg"),2), rep(("Otway"),2)),4),
         Treatment = rep(c("Unbaited", "Baited"),8))
# fix names 'cos of the 2.5 and 97.5 
names(det_preds) <- c("Predicted", "SD", "Lower", "Upper", "region_bait01", "XGROUPNAME", "Species", "Region", "Treatment") 

# reorder species
det_preds$Species <- fct_relevel(det_preds$Species, "Red fox", "Feral cat", "Southern brown bandicoot", "Long-nosed potoroo")

# plot
plot_daily_det <- ggplot(data=det_preds, aes(x=Region, y=Predicted, colour = Treatment)) +
  geom_pointrange(aes(ymin=Lower, ymax=Upper), position = position_dodge(width = 0.25), size = 0.3) +
  facet_wrap(~Species, scales = "free", nrow = 1) + 
  scale_color_manual(values = c("tomato", "dodgerblue"), labels = c("baited", "unbaited")) +   
  labs(subtitle = "a",x = "", y = "Daily Pr(detection)") +
  theme(plot.title = element_text(size=11),
        axis.title = element_text(size = 10),
        legend.position = "none",
        plot.margin = unit(c(0,0,0,0), "in"))
plot_daily_det

# plot / save
#png("figs/daily_detectability.png", width = 7, height = 7, res = 600, units = "in")
#plot_daily_det
#dev.off()


# PLOT: OCCUPANCY - FOX CONTROL ---------------------------------------------------------------

# now predict occupancy
fox_occ_pred <- predict(fit_occu_fox, newdata = df, submodel = "state", re.form = NA) %>% bind_cols(df)
cat_occ_pred <- predict(fit_occu_cat, newdata = df, submodel = "state", re.form = NA) %>% bind_cols(df)
sbb_occ_pred <- predict(fit_occu_sbb, newdata = df, submodel = "state", re.form = NA) %>% bind_cols(df)
lnp_occ_pred <- predict(fit_occu_lnp, newdata = df, submodel = "state", re.form = NA) %>% bind_cols(df)

# combine
occ_preds <- bind_rows(fox_occ_pred, cat_occ_pred, sbb_occ_pred, lnp_occ_pred) %>%
  mutate(Species = c(rep(("Red fox"),4), rep(("Feral cat"),4), rep(("Southern brown bandicoot"),4),  rep(("Long-nosed potoroo"),4)),
         Region = rep(c(rep(("Glenelg"),2), rep(("Otway"),2)),4),
         Treatment = rep(c("Unbaited", "Baited"),8)
         )
# fix names 'cos of the 2.5 and 97.5 
names(occ_preds) <- c("Predicted", "SD", "Lower", "Upper", "region_bait01", "XGROUPNAME", "Species", "Region", "Treatment") 

# reorder species
occ_preds$Species <- fct_relevel(occ_preds$Species, "Red fox", "Feral cat", "Southern brown bandicoot", "Long-nosed potoroo")

# plot
plot_occ <- ggplot(data=occ_preds, aes(x=Region, y=Predicted, colour = Treatment)) +
  geom_pointrange(aes(ymin=Lower, ymax=Upper), position = position_dodge(width = 0.25), size = 0.3) +
  facet_wrap(~Species, scales = "free", nrow = 1) + 
  scale_color_manual(values = c("tomato", "dodgerblue"), labels = c("baited", "unbaited")) +   
  labs(title = "",x = "Region", y = "Pr(occupancy)") +
  theme(plot.title = element_text(size=11),
        axis.title = element_text(size = 10),
        legend.position = "bottom",
        plot.margin = unit(c(0,0,0,0), "in"))

plot_occ



# PLOT: VEGETATION --------------------------------------------------------------

# make a dataframe to predict into
df <- expand.grid(region_bait01 = "otways_unbaited", XGROUPNAME = unique(records$XGROUPNAME))
# predict into
fox_veg_pred <- predict(fit_occu_fox, newdata = df, submodel = "state", re.form = NA) %>% bind_cols(df)
cat_veg_pred <- predict(fit_occu_cat, newdata = df, submodel = "state", re.form = NA) %>% bind_cols(df)
sbb_veg_pred <- predict(fit_occu_sbb, newdata = df, submodel = "state", re.form = NA) %>% bind_cols(df)
lnp_veg_pred <- predict(fit_occu_lnp, newdata = df, submodel = "state", re.form = NA) %>% bind_cols(df)

# combine
veg_preds <- bind_rows(fox_veg_pred, cat_veg_pred, sbb_veg_pred, lnp_veg_pred) %>%
  mutate(Species = c(rep(("Red fox"),7), rep(("Feral cat"),7), rep(("Southern brown bandicoot"),7),  rep(("Long-nosed potoroo"),7))
  )
# fix names 'cos of the 2.5 and 97.5 
names(veg_preds) <- c("Predicted", "SD", "Lower", "Upper", "region_bait01", "XGROUPNAME", "Species") 

# reorder species
veg_preds$Species <- fct_relevel(veg_preds$Species, "Red fox", "Feral cat", "Southern brown bandicoot", "Long-nosed potoroo")

# plot
plot_veg <- ggplot(data=veg_preds, aes(x=XGROUPNAME, y=Predicted)) +
  geom_pointrange(aes(ymin=Lower, ymax=Upper), size = 0.3) +
  facet_wrap(~Species, scales = "free", ncol = 1) + 
  labs(title = "",x = "Vegetation type", y = "Pr(occupancy)") +
  theme(plot.title = element_text(size=11),
        axis.title = element_text(size = 10),
        legend.position = "bottom",
        plot.margin = unit(c(0,0,0,0), "in"))
plot_veg

# plot / save
#png("figs/occupancy_vegetation.png", width = 8, height = 10, res = 600, units = "in")
#plot_veg
#dev.off()

# plot / save
png("figs/occ_det_fox_control.png", width = 8, height = 10.5, res = 600, units = "in")
plot_daily_det / plot_occ / plot_veg +
  plot_annotation(tag_levels = "a") +
  theme(plot.margin = margin(0, 0, 0, 0)) + 
  plot_layout(heights = unit(c(1, 1, 5), c('in', 'in', 'in')))
dev.off()

# END