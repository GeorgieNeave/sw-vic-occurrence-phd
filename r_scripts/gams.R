
rm(list = ls())
options(scipen = 999) 

## Load packages
library(mgcv)
# tools
library(dplyr)
library(forcats)
# model checking
library(DHARMa)


# LOAD DATA ------------------------------------------------------------

records <- read.csv("raw_data/raw_data_gams_pa.csv")

## VARIABLES
# region - Glenelg or Otway Ranges
# block - forest block. Glenelg = distinct forests, Otways = focal areas (Otway Ark data) + north / south (matts cat grids / zoi's surveys)
# data_source - dataset owner: Matt Rees (phd surveys), Zoi Banikos (masters surveys), Glenelg Ark camera-trap dataset 2013-19, Otway Ark dataset 2016-18
# station - unique ID for camera-trap site
# longitude - coordiantes for location of cam-trap site
# latitude - coordiantes for location of cam-trap site
# year - year of survey
# station_year - unique deployment: camera trap site x year 
# date_start - start date of cam-trap survey (yyyy-mm-dd)
# date_end - end date of camera-trap survey (yyyy-mm-dd)
# survey_duration - total survey duration (days) 
# fox - red fox: detected at least once (1) or not (0)
# cat - feral cat: detected at least once (1) or not (0)
# bandicoot_sb - bandicoot southern brown: detected at least once (1) or not (0)
# potoroo_ln - long-nosed potoroo: detected at least once (1) or not (0)
# XGROUPNAME - vegetation type: Ecological Vegetation Class group (DEWLP)
# EFG_NAME - vegetation type: Ecological Fire Group (https://www.ffm.vic.gov.au/__data/assets/pdf_file/0008/21113/Report-84-REDUCED-SIZE-Growth-Stages-and-Tolerable-Fire-Intervals-For-Victorias-Native-Vegetation-Data-Se.pdf)
# treatment - non-treatment / treatment (control / impact) sites in regards to fox control (i.e., blocks without fox control but were eventually subject to fox-baiting considered 'treatment') 
# foxbaits - number of 1080 fox-off poison-baits within a 2.3 km radius around each cam-trap (average fox dmax in Hradsky et al 2017 Sci Reports)
# elevation - metres above sea level based of Vicmap digital elevation model layer: 10m resolution (https://www.land.vic.gov.au/maps-and-spatial/spatial-data/vicmap-catalogue/vicmap-elevation)
# ruggedness - topographic ruggedness index made using QGIS plugin based off above^ elevation layer: median value within 30 m radius of cam-trap 
# wetness - topographic wetness index: median value within 30 m radius of cam-trap  (https://data.gov.au/dataset/ds-dap-csiro%3A5588/details?q=)
# rain_diff_percent_06months - deviation in rainfall (%) from the long-term average (closest weather station) in the last 6 months (from start of survey duration) 
# rain_diff_percent_12months - deviation in rainfall (%) from the long-term average (closest weather station) in the last 12 months (from start of survey duration) 
# rain_diff_percent_18months - deviation in rainfall (%) from the long-term average (closest weather station) in the last 18 months (from start of survey duration) 
# rain_diff_percent_24months - deviation in rainfall (%) from the long-term average (closest weather station) in the last 24 months (from start of survey duration) 
# tsf - time since last fire (year)
# n_fires - number of fires since 1939 bushfires
# dist_nnv - distance to nearest area of non-native vegatation (greater than 30 ha): derived from inverted extent of DELWP native vegetation extent layer


# ADJUST DATA ------------------------------------------------------------

## Drop camera-traps from dataset
# drop zoi's surveys - too much spatial autocorrelation (50 x 50 m grids)
records <- filter(records, data_source != "zoi")
# exclude stations left out for less than 10 days - can't trust em
records <- filter(records, survey_duration >= 10)
# only had three deployments of "Riverine Grassy Woodlands or Forests" vegetation type --> just remove em
records <- filter(records, XGROUPNAME != "Riverine Grassy Woodlands or Forests")

## Adjust vegetation type
# as we only had 20 "Rainforests" unique survey sites, all of which interspersed in "Wet Forests" in the Otways, reclassify as "Wet Forests". 
records$XGROUPNAME <- if_else(records$XGROUPNAME == "Rainforests", "Wet or Damp Forests", as.character(records$XGROUPNAME))
# abbreviate EVC group names for plotting
records$XGROUPNAME <- if_else(records$XGROUPNAME == "Wet or Damp Forests", "Wet Forests", records$XGROUPNAME)
records$XGROUPNAME <- if_else(records$XGROUPNAME == "Riparian Scrubs or Swampy Scrubs and Woodlands", "Swampy Scrubs", records$XGROUPNAME)
records$XGROUPNAME <- substr(records$XGROUPNAME, 1, nchar(records$XGROUPNAME) - 1)

## Adjust site random effect
## we need a random intercet for site to account for repeat sampling
## but 425 cam-traps in the Glenelg region deployed by Matt were only surveyed once (4 dense grids (~500 m spacing) in 4 forest blocks) - we can save 421 model parameters by renaming these with just the forest block
records$station <- if_else(records$data_source == "matt" & records$region == "glenelg", as.character(records$block), as.character(records$station))


## Check collinearity
# check correlation between numeric explanatory variables
cor <- as.data.frame(cor(records[,20:length(records)]))
# high correlation (-0.76) between time since fire (tsf) and number of historic fires (n_fires) --> keep in for now (one variable may be selected out), but pay attention to concurvity if not. 

## Transform variable classes
records <- transform(records,
                        region = factor(region, ordered = FALSE),                    
                        block = factor(block, ordered = FALSE),                      
                        data_source = factor(data_source, ordered = FALSE),                 
                        station = factor(station, ordered = FALSE),                     
                        station_year = factor(station_year, ordered = FALSE),             
                        survey_duration = as.integer(survey_duration),       
                        fox = as.integer(fox),                        
                        cat = as.integer(cat),                        
                        bandicoot_sb = as.integer(bandicoot_sb),              
                        potoroo_ln = as.integer(potoroo_ln),                 
                        XGROUPNAME = factor(XGROUPNAME, ordered = FALSE),              
                        foxbaits = round(as.numeric(foxbaits / 16.61902513749), 2), # convert to baits per km2              
                        elevation = as.integer(elevation),                
                        ruggedness = round(scale(as.numeric(ruggedness)), 2),                 
                        wetness = round(scale(as.numeric(wetness)), 2),                    
                        rain_diff_percent_06months = as.integer(rain_diff_percent_06months),  
                        rain_diff_percent_12months = as.integer(rain_diff_percent_12months),  
                        rain_diff_percent_18months = as.integer(rain_diff_percent_18months),  
                        rain_diff_percent_24months= as.integer(rain_diff_percent_24months),  
                        tsf = as.integer(tsf),                         
                        n_fires = as.integer(n_fires),                    
                        dist_nnv = as.integer(dist_nnv)
)



# FIT MODELS --------------------------------------------------------------

# Function to run the same model for each species / rainfall months  (we expect all variables to impact each species, and aim to contrast responses).
gam_function <- function(species, data, rainfall_months) {
  # rename column to species of interest
  data_renamed <- rename(data, species = species, rainfall_months = rainfall_months)
  # run GAM
  gam_model <- bam(species ~ s(tsf, bs = "tp", k = 7) +       
                             s(tsf, XGROUPNAME, bs = "fs", xt = list(bs = "tp"), k = 7) +
                           # s(n_fires, bs = "tp", k = 4) +   # removed due to concurvity issues
                            s(ruggedness, k = 4) +           # removed due to concurvity issues
                             s(wetness, k = 4) + 
                             s(foxbaits, region, bs = "fs", xt = list(bs = "tp"), k = 4) +
                             s(elevation, k = 4) + 
                             s(dist_nnv, bs = "tp", k = 4) +   
                             s(rainfall_months, region, bs = "fs", xt = list(bs = "tp"), k = 4) + 
                             s(station, bs = "re") +  
                             offset(log(survey_duration)), 
                   data = data_renamed, family = binomial, nthreads = 3, discrete = TRUE, select = TRUE)
  # save model
  return(gam_model)
}

## Fit GAM for each species x rainfall period
# fox
gam_fox_06 <- gam_function(species = "fox", data = records, rainfall_months = "rain_diff_percent_06months")
gam_fox_12 <- gam_function(species = "fox", data = records, rainfall_months = "rain_diff_percent_12months")
gam_fox_18 <- gam_function(species = "fox", data = records, rainfall_months = "rain_diff_percent_18months")
gam_fox_24 <- gam_function(species = "fox", data = records, rainfall_months = "rain_diff_percent_24months")
aic_fox <- AIC(gam_fox_06, gam_fox_12, gam_fox_18, gam_fox_24)

# cat
gam_cat_06 <- gam_function(species = "cat", data = records, rainfall_months = "rain_diff_percent_06months")
gam_cat_12 <- gam_function(species = "cat", data = records, rainfall_months = "rain_diff_percent_12months")
gam_cat_18 <- gam_function(species = "cat", data = records, rainfall_months = "rain_diff_percent_18months")
gam_cat_24 <- gam_function(species = "cat", data = records, rainfall_months = "rain_diff_percent_24months")
aic_cat <- AIC(gam_cat_06, gam_cat_12, gam_cat_18, gam_cat_24)

# southern brown bandicoot
gam_sbb_06 <- gam_function(species = "bandicoot_sb", data = records, rainfall_months = "rain_diff_percent_06months")
gam_sbb_12 <- gam_function(species = "bandicoot_sb", data = records, rainfall_months = "rain_diff_percent_12months")
gam_sbb_18 <- gam_function(species = "bandicoot_sb", data = records, rainfall_months = "rain_diff_percent_18months")
gam_sbb_24 <- gam_function(species = "bandicoot_sb", data = records, rainfall_months = "rain_diff_percent_24months")
aic_sbb <- AIC(gam_sbb_06, gam_sbb_12, gam_sbb_18, gam_sbb_24)

# long-nosed potoroo
gam_lnp_06 <- gam_function(species = "potoroo_ln", data = records, rainfall_months = "rain_diff_percent_06months")
gam_lnp_12 <- gam_function(species = "potoroo_ln", data = records, rainfall_months = "rain_diff_percent_12months")
gam_lnp_18 <- gam_function(species = "potoroo_ln", data = records, rainfall_months = "rain_diff_percent_18months")
gam_lnp_24 <- gam_function(species = "potoroo_ln", data = records, rainfall_months = "rain_diff_percent_24months")
aic_lnp <- AIC(gam_lnp_06, gam_lnp_12, gam_lnp_18, gam_lnp_24)


## Create table of AIC values
aic_df <- bind_rows(aic_fox, aic_cat, aic_sbb, aic_lnp) %>%
  mutate(species = c(rep("fox", 4), rep("cat", 4), rep("SBB", 4), rep("LNP", 4)),
         month = rep(c("6", "12", "18", "24"), 4))
# make a dAIC col
aic_df <- aic_df %>%
  group_by(species) %>%
  arrange(AIC, .by_group = TRUE) %>%
  mutate(dAIC = AIC - first(AIC)) %>%
  relocate(species, month, df, AIC, dAIC)
# save as csv
write.csv(aic_df, "derived_data/rain_occ_aic_table.csv")
## -> use 6 month rainfall for all species


## Check top models
summary(gam_fox_06)
summary(gam_cat_06)
summary(gam_sbb_06)
summary(gam_lnp_12)

# plot smooths
plot(gam_fox_06, seWithMean = TRUE, shade = TRUE, pages = 1, scheme = 2, scale = 0, rug = FALSE)
plot(gam_cat_06, seWithMean = TRUE, shade = TRUE, pages = 1, scheme = 2, scale = 0, rug = FALSE)
plot(gam_sbb_06, seWithMean = TRUE, shade = TRUE, pages = 1, scheme = 2, scale = 0, rug = FALSE)
plot(gam_lnp_12, seWithMean = TRUE, shade = TRUE, pages = 1, scheme = 2, scale = 0, rug = FALSE)

# check concurvity
concurvity(gam_fox_06, full = FALSE)
concurvity(gam_cat_06, full = FALSE)
concurvity(gam_sbb_06, full = FALSE)
concurvity(gam_lnp_12, full = FALSE)


# CHECK MODEL FITS AGAINST A NULL MODEL -----------------------------------
# how much are the random effects for site contributing to model fits?
# Function to run the null model 
gam_function_null <- function(species, data) {
  # rename column to species of interest
  data_renamed <- rename(data, species = species)
  # run GAM
  gam_model <- bam(species ~ s(station, bs = "re") +  
                     offset(log(survey_duration)), 
                   data = data_renamed, family = binomial, nthreads = 3, discrete = TRUE, select = TRUE)
  # save model
  return(gam_model)
}

# fit null
gam_fox_null <- gam_function_null(species = "fox", data = records)
gam_cat_null <- gam_function_null(species = "cat", data = records)
gam_sbb_null <- gam_function_null(species = "bandicoot_sb", data = records)
gam_lnp_null <- gam_function_null(species = "potoroo_ln", data = records)


# Model summary table -----------------------------------------------------
extract_fits <- function(model, species, model_name){
  x <- summary(model)
  df <- expand.grid(species = species,
                    model = model_name,
                    edf = sum(x$edf),
                    dev.expl = x$dev.expl,
                    r.sq = x$r.sq)
  return(df)
}

summaries <- bind_rows(extract_fits(gam_fox_null, "fox", "null"),
                       extract_fits(gam_cat_null, "cat", "null"),
                       extract_fits(gam_sbb_null, "sbb", "null"),
                       extract_fits(gam_lnp_null, "lnp", "null"),
                       extract_fits(gam_fox_06, "fox", "full"),
                       extract_fits(gam_cat_06, "cat", "full"),
                       extract_fits(gam_sbb_06, "sbb", "full"),
                       extract_fits(gam_lnp_12, "lnp", "full")) %>%
  arrange(species)

# add AIC scores
aic_vals <- AIC(gam_fox_null, gam_fox_06, gam_cat_null, gam_cat_06, gam_sbb_null, gam_sbb_06, gam_lnp_null, gam_lnp_12)
summaries <- cbind(summaries, aic_vals[2])

# save
write.csv(summaries, "derived_data/model_summaries.csv")


# --> plot scripts
# END