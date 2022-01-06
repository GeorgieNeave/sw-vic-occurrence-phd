library(ggplot2)
library(patchwork)

# Set the default theme for ggplot objects 
theme_set(theme_bw())
theme_update(panel.grid = element_blank())


# PLOT: TSF X VEG --------------------------------------------------------
gam_plot_tsf <- function(gam_model, data, title){
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
                  fitted = ilink(fit),
                  upper  = ilink(fit + (1.96 * se.fit)),
                  lower  = ilink(fit - (1.96 * se.fit)))
  # plot 
  plot_veg <- ggplot(data=df, aes(x=tsf, y=fitted, group=XGROUPNAME)) +
    geom_line(aes(y=fitted, x=tsf), lwd = 1.0) +
    facet_wrap(~XGROUPNAME, nrow = 1) +
    geom_ribbon(aes(ymin=lower, ymax=upper), alpha=0.2) +
    geom_rug(sides="t", data = filter(data, region == "glenelg"), aes(x = tsf), length = unit(0.05, "npc"), size = 0.2, col = "#440154FF", alpha = 0.8, inherit.aes = FALSE) +
    geom_rug(sides="b", data = filter(data, region == "otways"), aes(x = tsf),  length = unit(0.05, "npc"), size = 0.2, col = "#20A387FF", alpha = 0.8, inherit.aes = FALSE) +
    scale_color_manual(values = c("Glenelg" = "#440154FF", "Otway" = "#20A387FF"), name = "Region") +   
    scale_fill_manual(values = c("Glenelg" = "#440154FF", "Otway" = "#20A387FF"), name = "Region") +   
    labs(title = title, x = "Time since fire (years)", y = "Pr(occurrence)") + 
    theme(plot.title = element_text(size=11),
          axis.title = element_text(size = 10))
  plot_veg 
}

p1 <- gam_plot_tsf(gam_fox_06, records, expression(paste("Red fox ")))
p2 <- gam_plot_tsf(gam_cat_06, records, expression(paste("Feral cat ")))
p3 <- gam_plot_tsf(gam_sbb_06, records, expression(paste("Southern brown bandicoot ")))
p4 <- gam_plot_tsf(gam_lnp_12, records, expression(paste("Long-nosed potoroo ")))
tsf_plot <- p1 / p2 / p3 / p4 + plot_annotation(tag_levels = "a") 
tsf_plot

# plot / save
png("figs/tsf.png", width = 10, height = 9, res = 600, units = "in")
tsf_plot
dev.off()


# RAW DATA: TSF X VEG  -------------------------------------

# rename treatment variable
records$treatment <- ifelse(records$treatment == "treatment", "Baited", "Unbaited")

## Glenelg ark dataset
plot <- ggplot(filter(records, data_source == "glenelg_ark"), aes(x=XGROUPNAME, y=tsf, fill = XGROUPNAME)) + 
  geom_violin() + 
  facet_wrap(~treatment, nrow =1) +    
  scale_fill_discrete(drop=FALSE) + # retain unused levels (wet forest...)
  scale_x_discrete(drop=FALSE) +    # retain unused levels (wet forest...)
  labs(subtitle = "Glenelg Ark camera-trapping dataset (2013 - 2019)", x = "Vegetation type", y = "Time since fire (years)") + 
  theme(plot.title = element_text(size=11),
        axis.title = element_text(size = 10),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.position = "none",
        legend.title = element_blank())

## All datasets
plot_all <- ggplot(records, aes(x=XGROUPNAME, y=tsf, fill = XGROUPNAME)) + 
  geom_violin() + 
  facet_wrap(~treatment, nrow =1) +    
  labs(subtitle = "Entire dataset dataset (2013 - 2019)", x = "Vegetation type", y = "Time since fire (years)") + 
  theme(plot.title = element_text(size=11),
        axis.title = element_text(size = 10),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.position = "bottom",
        legend.title = element_blank())

# save
png("figs/raw_data_tsf_veg.png", width = 8, height = 10, res = 600, units = "in")
plot / plot_all
dev.off()


# END