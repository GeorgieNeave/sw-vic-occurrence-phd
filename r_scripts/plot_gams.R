
# SET-UP ------------------------------------------------------------------

library(ggplot2)
library(patchwork)
library(ggpubr)

# Set the default theme for ggplot objects 
theme_set(theme_bw())
theme_update(panel.grid = element_blank())

## remove outliers from plots
records$dist_nnv <- ifelse(records$dist_nnv > 5000, NA, records$dist_nnv)
records$ruggedness <- ifelse(records$ruggedness > 4, NA, records$ruggedness)
summary(records$ruggedness)

# LEGENDS -----------------------------------------------------------------
 # make plots to extract legend
gam_model <- gam_fox_06
data <- records

## VEG TYPE LEGEND
  # get model smooths
  x = sapply(gam_model$smooth, "[[",  "label")
  # make a new dataframe to predict into
  df <- expand.grid(XGROUPNAME = levels(data$XGROUPNAME),
                    survey_duration = 60, 
                    tsf = as.integer(seq(min(data$tsf), max(data$tsf), by = 1)))
  # adjust tsf min max based on actual survey range
  # first get real range for each veg group
  veg_tsf_range <- records %>%
    dplyr::group_by(XGROUPNAME) %>%
    mutate(tsf_min = min(tsf),
           tsf_max = max(tsf)) %>%
    select(XGROUPNAME, tsf_min, tsf_max) %>%
    distinct()
  # add to new dataframe and filter to correct range
  df <- left_join(df, veg_tsf_range) %>%
    filter(tsf>=tsf_min & tsf<=tsf_max)
  # predict model estimates (and uncertainty) into this dataframe, excluding the impact of variables other than hour:
  df <- cbind(df, predict.gam(gam_model, newdata = df, se.fit = TRUE, type = "link", newdata.guaranteed = TRUE, exclude = c(x[3:length(x)])))
  # get predictions on response scale
  ilink <- family(gam_model)$linkinv
  df <- transform(df,
                  upper  = ilink(fit + (1.96 * se.fit)),
                  lower  = ilink(fit - (1.96 * se.fit)),
                  fit = ilink(fit))
  # rename plots
  names(df)[1] <- "Vegetation_type"
  # plot
  plot_veg <- ggplot(data=df, aes(x=tsf, y=fit, colour = Vegetation_type, group = Vegetation_type)) +
    geom_line(aes(y=fit, x=tsf), lwd = 0.75) +
    geom_ribbon(aes(ymin=lower, ymax=upper, fill = Vegetation_type), lty = 0, alpha=0.075) +
    geom_rug(sides="t", data = filter(data, region == "glenelg"), aes(x = tsf), length = unit(0.04, "npc"), size = 0.2, col = "#440154FF", alpha = 0.8, inherit.aes = FALSE) +
    geom_rug(sides="b", data = filter(data, region == "otways"), aes(x = tsf),  length = unit(0.04, "npc"), size = 0.2, col = "#20A387FF", alpha = 0.8, inherit.aes = FALSE) +
    labs(title = "Time since fire: vegetation type", x = "Years", y = "") + 
    theme(plot.title = element_text(size=11),
          axis.title = element_text(size = 10),
          legend.position = "right",
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          plot.margin = unit(c(0,0,0,0), "in"))
# extract legend (via ggpubr)
legend_veg <- as_ggplot(get_legend(plot_veg, position = NULL))
  
  
## REGION LEGEND
  # make a new dataframe to predict into
  df <- expand.grid(survey_duration = 60,
                    region = c("glenelg", "otways"),
                    foxbaits = seq(min(data$foxbaits), max(data$foxbaits), by = 0.05))
  # predict model estimates (and uncertainty) into this dataframe, excluding the impact of other variables
  df <- cbind(df, predict.gam(gam_model, newdata = df, se.fit = TRUE, type = "link", newdata.guaranteed = TRUE, exclude = c(x[1:4], x[6:length(x)])))
  # get predictions on response scale
  df <- transform(df,
                  upper  = ilink(fit + (1.96 * se.fit)),
                  lower  = ilink(fit - (1.96 * se.fit)),
                  fit = ilink(fit))
  # plot
plot_fc <-  ggplot(data=df, aes(x=foxbaits, y=fit, col=region)) +
    geom_ribbon(aes(ymin=lower, ymax=upper, fill=region), alpha=0.2, linetype=0) +
    geom_line(aes(y=fit, x=foxbaits), lwd = 1.0) +
    scale_fill_manual(values = c("#440154FF", "#20A387FF"), labels = c("Glenelg", "Otway"), name = "Region") +
    scale_color_manual(values = c("#440154FF", "#20A387FF"), labels = c("Glenelg", "Otway"), name = "Region") +   
    geom_rug(sides="t", data = filter(data, region == "glenelg"), aes(x = foxbaits), length = unit(0.04, "npc"), size = 0.2, col = "#440154FF", alpha = 0.8, inherit.aes = FALSE) +
    geom_rug(sides="b", data = filter(data, region == "otways"), aes(x = foxbaits),  length = unit(0.04, "npc"), size = 0.2, col = "#20A387FF", alpha = 0.8, inherit.aes = FALSE) +
    labs(title = "Fox control", x = bquote('1080 poison ' (baits ~km^-2)), y = "Pr(occurrence)") + 
    theme(plot.title = element_text(size=11),
          axis.title = element_text(size = 10),
          legend.position = "bottom",
          plot.margin = unit(c(0,0,0,0), "in"))

# extract legend / save
legend_reg <- as_ggplot(get_legend(plot_fc, position = NULL))





# 1) TSF GS ---------------------------------------------------------------
func_tsf_gs <- function(gam_model, data, min, max){
  # get model smooths
  x = sapply(gam_model$smooth, "[[",  "label")
  # make a new dataframe to predict into
  df <- expand.grid(survey_duration = 60, 
                    tsf = as.integer(seq(min(data$tsf), max(data$tsf), by = 1)))
  # predict model estimates :
  df <- cbind(df, predict.gam(gam_model, newdata = df, se.fit = TRUE, type = "link", newdata.guaranteed = TRUE, exclude = c(x[2:length(x)])))
  # get predictions on response scale
  ilink <- family(gam_model)$linkinv
  df <- transform(df,
                  upper  = ilink(fit + (1.96 * se.fit)),
                  lower  = ilink(fit - (1.96 * se.fit)),
                  fit = ilink(fit))
  # plot
  ggplot(data=df, aes(x=tsf, y=fit)) +
    geom_line(aes(y=fit, x=tsf), lwd = 1.0) +
    ylim(min, max) +
    geom_ribbon(aes(ymin=lower, ymax=upper), alpha=0.2) +
    geom_rug(sides="t", data = filter(data, region == "glenelg"), aes(x = tsf), length = unit(0.04, "npc"), size = 0.2, col = "#440154FF", alpha = 0.8, inherit.aes = FALSE) +
    geom_rug(sides="b", data = filter(data, region == "otways"), aes(x = tsf),  length = unit(0.04, "npc"), size = 0.2, col = "#20A387FF", alpha = 0.8, inherit.aes = FALSE) +
    labs(title = "Time since fire: average", x = "Years", y = "Pr(occurrence)") + 
    theme(plot.title = element_text(size=11),
          axis.title = element_text(size = 10),
          plot.margin = unit(c(0,0,0,0), "in"))
}


# 2) TSF X VEG -------------------------------------------------------------
func_tsf_veg <- function(gam_model, data, min, max){
  # get model smooths
  x = sapply(gam_model$smooth, "[[",  "label")
  # make a new dataframe to predict into
  df <- expand.grid(XGROUPNAME = levels(data$XGROUPNAME),
                    survey_duration = 60, 
                    tsf = as.integer(seq(min(data$tsf), max(data$tsf), by = 1)))
  # adjust tsf min max based on actual survey range
  # first get real range for each veg group
  veg_tsf_range <- records %>%
    dplyr::group_by(XGROUPNAME) %>%
    mutate(tsf_min = min(tsf),
           tsf_max = max(tsf)) %>%
    select(XGROUPNAME, tsf_min, tsf_max) %>%
    distinct()
  # add to new dataframe and filter to correct range
  df <- left_join(df, veg_tsf_range) %>%
    filter(tsf>=tsf_min & tsf<=tsf_max)
  # predict model estimates (and uncertainty) into this dataframe, excluding the impact of variables other than hour:
  df <- cbind(df, predict.gam(gam_model, newdata = df, se.fit = TRUE, type = "link", newdata.guaranteed = TRUE, exclude = c(x[3:length(x)])))
  # get predictions on response scale
  ilink <- family(gam_model)$linkinv
  df <- transform(df,
                  upper  = ilink(fit + (1.96 * se.fit)),
                  lower  = ilink(fit - (1.96 * se.fit)),
                  fit = ilink(fit))
  # rename plots
  names(df)[1] <- "Vegetation_type"
  # plot
 plot_veg <- ggplot(data=df, aes(x=tsf, y=fit, colour = Vegetation_type, group = Vegetation_type)) +
    ylim(min, max) +
    geom_line(aes(y=fit, x=tsf), lwd = 0.75) +
    geom_ribbon(aes(ymin=lower, ymax=upper, fill = Vegetation_type), lty = 0, alpha=0.075) +
    geom_rug(sides="t", data = filter(data, region == "glenelg"), aes(x = tsf), length = unit(0.04, "npc"), size = 0.2, col = "#440154FF", alpha = 0.8, inherit.aes = FALSE) +
    geom_rug(sides="b", data = filter(data, region == "otways"), aes(x = tsf),  length = unit(0.04, "npc"), size = 0.2, col = "#20A387FF", alpha = 0.8, inherit.aes = FALSE) +
    labs(title = "Time since fire: vegetation type", x = "Years", y = "") + 
    theme(plot.title = element_text(size=11),
          axis.title = element_text(size = 10),
          legend.position = "none",
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          plot.margin = unit(c(0,0,0,0), "in"))
}


# 3) TWI ------------------------------------------------------------
func_twi <- function(gam_model, data, min, max){
  # get model smooths
  x = sapply(gam_model$smooth, "[[",  "label")
  # make a new dataframe to predict into
  df <- expand.grid(survey_duration = 60, 
                    wetness = seq(min(data$wetness), max(data$wetness), by = 0.1))
  # predict model estimates :
  df <- cbind(df, predict.gam(gam_model, newdata = df, se.fit = TRUE, type = "link", newdata.guaranteed = TRUE, exclude = c(x[1:3], x[5:length(x)])))
  # get predictions on response scale
  ilink <- family(gam_model)$linkinv
  df <- transform(df,
                  upper  = ilink(fit + (1.96 * se.fit)),
                  lower  = ilink(fit - (1.96 * se.fit)),
                  fit = ilink(fit))
  # plot
  ggplot(data=df, aes(x=wetness, y=fit)) +
    ylim(min, max) +
    geom_line(aes(y=fit, x=wetness), lwd = 1.0) +
    geom_ribbon(aes(ymin=lower, ymax=upper), alpha=0.2) +
    geom_rug(sides="t", data = filter(data, region == "glenelg"), aes(x = wetness), length = unit(0.04, "npc"), size = 0.2, col = "#440154FF", alpha = 0.8, inherit.aes = FALSE) +
    geom_rug(sides="b", data = filter(data, region == "otways"), aes(x = wetness),  length = unit(0.04, "npc"), size = 0.2, col = "#20A387FF", alpha = 0.8, inherit.aes = FALSE) +
    labs(title = "Topographic wetness", x = "Index", y = "") + 
    theme(plot.title = element_text(size=11),
          axis.title = element_text(size = 10),
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          plot.margin = unit(c(0,0,0,0), "in"))
}


# 4) FOXBAITS -------------------------------------------------------------
func_foxbaits <- function(gam_model, data, min, max){
  # get model smooths
  x = sapply(gam_model$smooth, "[[",  "label")
  # make a new dataframe to predict into
  df <- expand.grid(survey_duration = 60,
                    rainfall_months = 0, 
                    region = c("glenelg", "otways"),
                    foxbaits = seq(min(data$foxbaits), max(data$foxbaits), by = 0.05))
  # predict model estimates (and uncertainty) into this dataframe, excluding the impact of other variables
  df <- cbind(df, predict.gam(gam_model, newdata = df, se.fit = TRUE, type = "link", newdata.guaranteed = TRUE, exclude = c(x[1:4], x[6:7], x[9:length(x)])))
  # get predictions on response scale
  ilink <- family(gam_model)$linkinv
  df <- transform(df,
                  upper  = ilink(fit + (1.96 * se.fit)),
                  lower  = ilink(fit - (1.96 * se.fit)),
                  fit = ilink(fit))
# plot
 plot_fc <-  ggplot(data=df, aes(x=foxbaits, y=fit, col=region)) +
    ylim(min, max) +
    geom_ribbon(aes(ymin=lower, ymax=upper, fill=region), alpha=0.2, linetype=0) +
    geom_line(aes(y=fit, x=foxbaits), lwd = 1.0) +
    scale_fill_manual(values = c("#440154FF", "#20A387FF"), labels = c("Glenelg", "Otway"), name = "Region") +
    scale_color_manual(values = c("#440154FF", "#20A387FF"), labels = c("Glenelg", "Otway"), name = "Region") +   
    geom_rug(sides="t", data = filter(data, region == "glenelg"), aes(x = foxbaits), length = unit(0.04, "npc"), size = 0.2, col = "#440154FF", alpha = 0.8, inherit.aes = FALSE) +
    geom_rug(sides="b", data = filter(data, region == "otways"), aes(x = foxbaits),  length = unit(0.04, "npc"), size = 0.2, col = "#20A387FF", alpha = 0.8, inherit.aes = FALSE) +
    labs(title = "Fox control", x = bquote('1080 poison ' (baits ~km^-2)), y = "Pr(occurrence)") + 
    theme(plot.title = element_text(size=11),
          axis.title = element_text(size = 10),
          legend.position = "none",
          plot.margin = unit(c(0,0,0,0), "in"))
}


# 5) ELEVATION ------------------------------------------------------------
func_elevation <- function(gam_model, data, min, max){
  # get model smooths
  x = sapply(gam_model$smooth, "[[",  "label")
  # make a new dataframe to predict into
  df <- expand.grid(survey_duration = 60, 
                    elevation = seq(min(data$elevation), max(data$elevation), by = 10))
  # predict model estimates :
  df <- cbind(df, predict.gam(gam_model, newdata = df, se.fit = TRUE, type = "link", newdata.guaranteed = TRUE, exclude = c(x[1:5], x[7:length(x)])))
  # get predictions on response scale
  ilink <- family(gam_model)$linkinv
  df <- transform(df,
                  upper  = ilink(fit + (1.96 * se.fit)),
                  lower  = ilink(fit - (1.96 * se.fit)),
                  fit = ilink(fit))
  # plot
  ggplot(data=df, aes(x=elevation, y=fit)) +
    ylim(min, max) +
    geom_line(aes(y=fit, x=elevation), lwd = 1.0) +
    geom_ribbon(aes(ymin=lower, ymax=upper), alpha=0.2) +
    geom_rug(sides="t", data = filter(data, region == "glenelg"), aes(x = elevation), length = unit(0.04, "npc"), size = 0.2, col = "#440154FF", alpha = 0.8, inherit.aes = FALSE) +
    geom_rug(sides="b", data = filter(data, region == "otways"), aes(x = elevation),  length = unit(0.04, "npc"), size = 0.2, col = "#20A387FF", alpha = 0.8, inherit.aes = FALSE) +
    labs(title = "Elevation", x = "Above sea level (metres)", y = "") + 
    theme(plot.title = element_text(size=11),
          axis.title = element_text(size = 10),
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          plot.margin = unit(c(0,0,0,0), "in"))
}



# 6) RUGGEDNESS ------------------------------------------------------------
func_ruggedness <- function(gam_model, data, min, max){
  # get model smooths
  x = sapply(gam_model$smooth, "[[",  "label")
  # make a new dataframe to predict into
  df <- expand.grid(survey_duration = 60, 
                    ruggedness = seq(-0.910, 4, by = 0.1))
  # predict model estimates :
  df <- cbind(df, predict.gam(gam_model, newdata = df, se.fit = TRUE, type = "link", newdata.guaranteed = TRUE, exclude = c(x[1:2], x[4:length(x)])))
  # get predictions on response scale
  ilink <- family(gam_model)$linkinv
  df <- transform(df,
                  upper  = ilink(fit + (1.96 * se.fit)),
                  lower  = ilink(fit - (1.96 * se.fit)),
                  fit = ilink(fit))
  # plot
  ggplot(data=df, aes(x=ruggedness, y=fit)) +
    ylim(min, max) +
    geom_line(aes(y=fit, x=ruggedness), lwd = 1.0) +
    geom_ribbon(aes(ymin=lower, ymax=upper), alpha=0.2) +
    geom_rug(sides="t", data = filter(data, region == "glenelg"), aes(x = ruggedness), length = unit(0.04, "npc"), size = 0.2, col = "#440154FF", alpha = 0.8, inherit.aes = FALSE) +
    geom_rug(sides="b", data = filter(data, region == "otways"), aes(x = ruggedness),  length = unit(0.04, "npc"), size = 0.2, col = "#20A387FF", alpha = 0.8, inherit.aes = FALSE) +
    labs(title = "Terrain ruggedness", x = "Index", y = "Pr(occurrence)") + 
    theme(plot.title = element_text(size=11),
          axis.title = element_text(size = 10),
          plot.margin = unit(c(0,0,0,0), "in"))
}


# 7) DIST NNV -------------------------------------------------------------
func_dist_nnv <- function(gam_model, data, min, max){
  # get model smooths
  x = sapply(gam_model$smooth, "[[",  "label")
  # make a new dataframe to predict into
  df <- expand.grid(survey_duration = 60, 
                    dist_nnv = seq(0, 5000, by = 100))
  # predict model estimates :
  df <- cbind(df, predict.gam(gam_model, newdata = df, se.fit = TRUE, type = "link", newdata.guaranteed = TRUE, exclude = c(x[1:6], x[8:length(x)])))
  # get predictions on response scale
  ilink <- family(gam_model)$linkinv
  df <- transform(df,
                  upper  = ilink(fit + (1.96 * se.fit)),
                  lower  = ilink(fit - (1.96 * se.fit)),
                  fit = ilink(fit))
  # convert to km
  df$dist_nnv <- df$dist_nnv / 1000
  # plot
  ggplot(data=df, aes(x=dist_nnv, y=fit)) +
    ylim(min, max) +
    geom_line(aes(y=fit, x=dist_nnv), lwd = 1.0) +
    geom_ribbon(aes(ymin=lower, ymax=upper), alpha=0.2) +
    geom_rug(sides="t", data = filter(data, region == "glenelg"), aes(x = dist_nnv/1000), length = unit(0.04, "npc"), size = 0.2, col = "#440154FF", alpha = 0.8, inherit.aes = FALSE) +
    geom_rug(sides="b", data = filter(data, region == "otways"), aes(x = dist_nnv/1000),  length = unit(0.04, "npc"), size = 0.2, col = "#20A387FF", alpha = 0.8, inherit.aes = FALSE) +
    labs(title = "Distance to forest edge", x = "Kilometres", y = "") + 
    theme(plot.title = element_text(size=11),
          axis.title = element_text(size = 10),
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          plot.margin = unit(c(0,0,0,0), "in"))
}


# 8) RECENT RAINFALL ------------------------------------------------------
func_rain <- function(gam_model, data, rainfall_months, month, min, max){
  # rename based on months
  data <- rename(data, rainfall_months = rainfall_months)
  # get model smooths
  x = sapply(gam_model$smooth, "[[",  "label")
  # make a new dataframe to predict into
  df <- expand.grid(foxbaits = mean(data$foxbaits), 
                    survey_duration = max(data$survey_duration), 
                    region = c("glenelg", "otways"),
                    rainfall_months = seq(min(data$rainfall_months), max(data$rainfall_months), by = 1))
  # adjust tsf min max based on actual survey range
  # first get real range for each veg group
  rain_range <- data %>%
    dplyr::group_by(region) %>%
    mutate(rainfall_months_min = min(rainfall_months),
           rainfall_months_max = max(rainfall_months)) %>%
    select(region, rainfall_months_min, rainfall_months_max) %>%
    distinct()
  # add to new dataframe and filter to correct range
  df <- left_join(df, rain_range) %>%
    filter(rainfall_months>=rainfall_months_min & rainfall_months<=rainfall_months_max)
  # predict model estimates (and uncertainty) into this dataframe, excluding the impact of variables other than hour:
  df <- cbind(df, predict.gam(gam_model, newdata = df, se.fit = TRUE, type = "link", newdata.guaranteed = TRUE, exclude = c(x[1:4], x[6:7], x[9])))
  # get predictions on response scale
  ilink <- family(gam_model)$linkinv
  df <- transform(df,
                  upper  = ilink(fit + (1.96 * se.fit)),
                  lower  = ilink(fit - (1.96 * se.fit)),
                  fit = ilink(fit))
# plot 
  ggplot(data=df, aes(x=rainfall_months, y=fit, col=region)) +
    ylim(min, max) +
    geom_ribbon(aes(ymin=lower, ymax=upper, fill=region), alpha=0.2, linetype=0) +
    geom_line(aes(y=fit, x=rainfall_months), lwd = 1.0) +
    scale_fill_manual(values = c("#440154FF", "#20A387FF"), labels = c("Glenelg", "Otway"), name = "Region") +
    scale_color_manual(values = c("#440154FF", "#20A387FF"), labels = c("Glenelg", "Otway"), name = "Region") +   
    geom_rug(sides="t", data = filter(data, region == "glenelg"), aes(x = rainfall_months), length = unit(0.04, "npc"), size = 0.2, col = "#440154FF", alpha = 0.8, inherit.aes = FALSE) +
    geom_rug(sides="b", data = filter(data, region == "otways"), aes(x = rainfall_months),  length = unit(0.04, "npc"), size = 0.2, col = "#20A387FF", alpha = 0.8, inherit.aes = FALSE) +
    geom_vline(xintercept = 0, colour = "black", size = 0.6, linetype="dotted") + 
    labs(title = "Recent rainfall", x = paste("Deviation (%)"), y = "") + 
    theme(plot.title = element_text(size=11),
        axis.title = element_text(size = 10),
        legend.position = "none",
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        plot.margin = unit(c(0,0,0,0), "in"))
}


# 9) REGION ------------------------------------------------------------
func_region <- function(gam_model, data, min, max){
  # get model smooths
  x = sapply(gam_model$smooth, "[[",  "label")
  # make a new dataframe to predict into
  df <- expand.grid(survey_duration = 60, 
                    foxbaits = 0,
                    rainfall_months = 0,
                    region = levels(data$region))
  # predict model estimates :
  df <- cbind(df, predict.gam(gam_model, newdata = df, se.fit = TRUE, type = "link", newdata.guaranteed = TRUE, exclude = c(x[1:4], x[5:7], x[9])))
  # get predictions on response scale
  ilink <- family(gam_model)$linkinv
  df <- transform(df,
                  upper  = ilink(fit + (1.96 * se.fit)),
                  lower  = ilink(fit - (1.96 * se.fit)),
                  fit = ilink(fit))
  # rename region
  df$region <- gsub("glenelg", "Glenelg", df$region)
  df$region <- gsub("otways", "Otway", df$region)
  # plot
  ggplot(data=df, aes(x=region, y=fit, col=region)) +
    ylim(min, max) +
    geom_point(size = 2) +
    geom_errorbar(aes(x=region, ymin=lower, ymax=upper), size = 0.5, width = 0.2) +
    scale_fill_manual(values = c("#440154FF", "#20A387FF"), labels = c("Glenelg", "Otway"), name = "Region") +
    scale_color_manual(values = c("#440154FF", "#20A387FF"), labels = c("Glenelg", "Otway"), name = "Region") +   
    labs(title = "Region", x = "", y = "Pr(occurrence)") + 
    theme(plot.title = element_text(size=11),
          axis.title = element_text(size = 10),
          legend.position = "none",
          plot.margin = unit(c(0,0,0,0), "in"))
}


# A) FOX PLOTS ------------------------------------------------------------
# run the functions
a <- func_region(gam_fox_06, records, 0, 0.78)
b <- func_foxbaits(gam_fox_06, records, 0, 0.78)
c <- func_rain(gam_fox_06, records, rainfall_months = "rain_diff_percent_06months", month = "24", 0, 0.78)
d <- func_elevation(gam_fox_06, records, 0, 0.78)
e <- func_ruggedness(gam_fox_06, records, 0, 0.78)
f <- func_twi(gam_fox_06, records, 0, 0.78)
g <- func_dist_nnv(gam_fox_06, records, 0, 0.78)
h <- func_tsf_gs(gam_fox_06, records, 0, 0.78)
i <- func_tsf_veg(gam_fox_06, records, 0, 0.78)

# assemble / save
png("figs/gams_fox.png", width = 8, height = 8.5, res = 600, units = "in")
( b | c | d ) / ( e | f | g) / ( h | i | legend_veg) + 
  plot_annotation(title = "Red fox", 
                  theme = theme(plot.title = element_text(size = 18)),  
                  tag_levels = list(c("a", "b", "c", "d", "e", "f", "g", "h", "")))  
dev.off()

# to combine region legend, type in the terminal (using imagemagick): 
# convert figs/gams_fox.png figs/gams_legend_reg.png  -append figs/gams_fox.png



# B) CAT PLOTS ------------------------------------------------------------
# run the functions
a <- func_region(gam_cat_06, records, 0, 0.65)
b <- func_foxbaits(gam_cat_06, records, 0, 0.65)
c <- func_rain(gam_cat_06, records, rainfall_months = "rain_diff_percent_06months", month = "24", 0, 0.65)
d <- func_elevation(gam_cat_06, records, 0, 0.65)
e <- func_ruggedness(gam_cat_06, records, 0, 0.65)
f <- func_twi(gam_cat_06, records, 0, 0.65)
g <- func_dist_nnv(gam_cat_06, records, 0, 0.65)
h <- func_tsf_gs(gam_cat_06, records, 0, 0.65)
i <- func_tsf_veg(gam_cat_06, records, 0, 0.65)
# assemble / save
png("figs/gams_cat.png", width = 8, height = 8.5, res = 600, units = "in")
( b | c | d ) / ( e | f | g) / ( h | i | legend_veg) + 
  plot_annotation(title = "Feral cat", 
                  theme = theme(plot.title = element_text(size = 18)),  
                  tag_levels = list(c("a", "b", "c", "d", "e", "f", "g", "h", "")))  
dev.off()

# to combine region legend, type in the terminal (using imagemagick): 
# convert figs/gams_cat.png figs/gams_legend_reg.png  -append figs/gams_cat.png


# C) SBB PLOTS ------------------------------------------------------------
# run the functions
a <- func_region(gam_sbb_06, records, 0, 0.2)
b <- func_foxbaits(gam_sbb_06, records, 0, 0.2)
c <- func_rain(gam_sbb_06, records, rainfall_months = "rain_diff_percent_06months", month = "24", 0, 0.2)
d <- func_elevation(gam_sbb_06, records, 0, 0.2)
e <- func_ruggedness(gam_sbb_06, records, 0, 0.2)
f <- func_twi(gam_sbb_06, records, 0, 0.2)
g <- func_dist_nnv(gam_sbb_06, records, 0, 0.2)
h <- func_tsf_gs(gam_sbb_06, records, 0, 0.2)
i <- func_tsf_veg(gam_sbb_06, records, 0, 0.2)

# assemble / save
png("figs/gams_sbb.png", width = 8, height = 8.5, res = 600, units = "in")
( b | c | d ) / ( e | f | g) / ( h | i | legend_veg) + 
  plot_annotation(title = "Southern brown bandicoot", 
                  theme = theme(plot.title = element_text(size = 18)),  
                  tag_levels = list(c("a", "b", "c", "d", "e", "f", "g", "h", "")))  
dev.off()

# to combine region legend, type in the terminal (using imagemagick): 
# convert figs/gams_sbb.png figs/gams_legend_reg.png  -append figs/gams_sbb.png


# D) LNP PLOTS ------------------------------------------------------------
# run the functions
a <- func_region(gam_lnp_12, records, 0, 0.7)
b <- func_foxbaits(gam_lnp_12, records, 0, 0.7)
c <- func_rain(gam_lnp_12, records, rainfall_months = "rain_diff_percent_12months", month = "24", 0, 0.7)
d <- func_elevation(gam_lnp_12, records, 0, 0.7)
e <- func_ruggedness(gam_lnp_12, records, 0, 0.7)
f <- func_twi(gam_lnp_12, records, 0, 0.7)
g <- func_dist_nnv(gam_lnp_12, records, 0, 0.7)
h <- func_tsf_gs(gam_lnp_12, records, 0, 0.7)
i <- func_tsf_veg(gam_lnp_12, records, 0, 0.7)

# assemble / save
png("figs/gams_lnp.png", width = 8, height = 8.5, res = 600, units = "in")
( b | c | d ) / ( e | f | g) / ( h | i | legend_veg) + 
  plot_annotation(title = "Long-nosed potoroo", 
                  theme = theme(plot.title = element_text(size = 18)),  
                  tag_levels = list(c("a", "b", "c", "d", "e", "f", "g", "h", "")))  
dev.off()

# to combine region legend, type in the terminal (using imagemagick): 
# convert figs/gams_lnp.png figs/gams_legend_reg.png  -append figs/gams_lnp.png


# END