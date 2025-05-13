# Author: Brandi Moore
# Date code compiled: 05-09-2025
# R Version: R 4.3.2
# Institution: The Center for Opioid Epidemiology and Policy (NYU Langone)

###########################################################
### Code to replicate primary analysis for:             ###
### Considerations for the epidemiological evaluation   ###
### of hyperlocal interventions: A case study of the    ###
### New York City overdose prevention centers           ### 
###########################################################

###########################################################
### CODE FOR SYNTHETIC CONTROL ANALYSIS WITH            ###
### 500M BUFFERS AROUND MOBILE/STATIONARY SSPs AS DONORS###
###########################################################

# load packages
library(sf)
library(measurements)
library(dplyr)
library(areal)
library(Synth)
library(doParallel)

######################################
########### CREATE BUFFERS ##########
######################################

# load in ssp points, opc points, and nyc census tract shp
ssp <- st_read("ssps.shp")
eharlem_opc <- st_read("harlemOPC.shp")
wheights_opc <- st_read("heightsOPC.shp")
nyc <- st_read("nyct.shp")

# spatial transformation to ensure matching crs with nyc ct shp
CRS.new <- st_crs("EPSG:2263")
eharlem_opc.2263 <- st_transform(eharlem_opc, CRS.new)
wheights_opc.2263 <- st_transform(wheights_opc, CRS.new)
ssp.2263 <- st_transform(ssp, CRS.new)

# check new crs
st_crs(eharlem_opc.2263) 
st_crs(wheights_opc.2263) 
st_crs(ssp.2263)

# generate 500 meter circular buffers

# convert meters to feet 
conv_unit(500,"m","ft")

# create buffers
eharlem_opc_buffer <- st_buffer(eharlem_opc.2263, 1640.42)
wheights_opc_buffer <- st_buffer(wheights_opc.2263, 1640.42)
ssp_buffer <- st_buffer(ssp.2263, 1640.42)

# merge buffers
buffers <- rbind(eharlem_opc_buffer, wheights_opc_buffer)
buffers <- rbind(buffers, ssp_buffer)
# create ID number for buffers
buffers$ID <- 1:nrow(buffers)

######################################
##### GET DAILY COUNTS OF ##########
##### DRUG-RELATED ARRESTS IN BUFFER 
######################################

# create initial blank df with correct number of rows to represent days of time period (01/01/2014 to 09/30/2023)
ID <- rep(1:28, each=3560)
day_numb <- rep(1:3560, times=28)
blank_df <- data.frame(ID, day_numb)

# merge blank df with buffers
buffer_blank <- merge(blank_df, buffers, by="ID")

# load in shp for drug arrests
drug_arrest.2263 <- st_read("drug_arrest_2023.sf.shp")

# get count of daily drug arrests within buffers
buffer_drugarrest <- st_intersection(x = buffers, y = drug_arrest.2263)
buffer_drugarrest <- buffer_drugarrest %>% add_count(day_numb, ID, name = 'buffer_daily_count')

# merge buffer repeat file and drug arrests buffer 
buffer_blank_drugarrests <- merge(buffer_blank, buffer_drugarrest, by=c("ID","day_numb", "Name", "Address","City","Zip","Full_Add", "State", "geo_method"), all=TRUE)

# make actual count of drug arrests per day per site 
arrest_count <- buffer_blank_drugarrests %>% 
  group_by(ID, day_numb, Name, Address, City, Zip, Full_Add, State, geo_method) %>% 
  summarise(tot_day_count = sum(!is.na(buffer_daily_count)))

# Merge actual counts per day per site with buffers
daily_arrest_count_buffer <- merge(buffer_blank, arrest_count, by=c("ID","day_numb", "Name", "Address","City","Zip","Full_Add", "State", "geo_method"))

# save data file if desired

######################################
## INTERPOLATE ACS VARS TO BUFFERS ###
######################################

# load in acs data
acs_2019 <- read.csv("ACS2019_5yr_CensusTract.csv")

# transform acs variables to numeric form
acs_2019[3:36] <- lapply(acs_2019[3:36], as.numeric)

# turn NAs -> 0
acs_2019[is.na(acs_2019)] <- 0

# merge acs with nyc on common variable 
acs_nyc <- merge(nyc, acs_2019, by="GEOID")

# interpolate acs to buffers
interpolated_acs_buffers <- aw_interpolate(buffers, tid = ID, source = acs_nyc, sid = GEOID, weight = "sum", output = "sf", 
                                           intensive = c("pct_male",
                                                         "pct_female", 
                                                         "median_household_income_in_past_12_months",
                                                         "pct_poverty",
                                                         "pct_unemployed",
                                                         "pct_public_assistance",
                                                         "median_age",
                                                         "pct_hs_ged_attainment",
                                                         "pct_ai_an",
                                                         "pct_asian",
                                                         "pct_black",
                                                         "pct_latino",
                                                         "pct_nh_pi",
                                                         "pct_white",
                                                         "pct_other",
                                                         "pct_two_or_more_race"),
                                           extensive = c("number_of_residents", 
                                                         "number_of_males", 
                                                         "number_of_females", 
                                                         "number_of_residents_with_income_below_poverty_level",
                                                         "number_of_residents_unemployed",
                                                         "number_of_residents_in_labor_force",
                                                         "number_of_households_with_public_assistance_income_in_past_12_months",
                                                         "number_of_households",
                                                         "number_of_residents_education_at_least_hsdegree_ged",
                                                         "number_of_residents_aged_25_older",
                                                         "number_of_white_residents",
                                                         "number_of_black_african_american_residents",
                                                         "number_of_american_indian_alaska_native_residents",
                                                         "number_of_asian_residents",
                                                         "number_of_native_hawaiian_pacific_islander_residents",
                                                         "number_of_other_race_residents",
                                                         "number_of_two_or_more_race_residents",
                                                         "number_of_hispanic_latino_residents"))

# prep for merging interpolated buffers with buffer-arrest data 

#replicate values of interpolated data
daily_acs <- interpolated_acs_buffers[rep(seq_len(nrow(interpolated_acs_buffers)), 3560),]

# add day_numb
daily_acs$day_numb <- ave(daily_acs$Name, daily_acs$Name, FUN=seq_along)

# ensure these vars are numeric
daily_acs$day_numb <- as.numeric(daily_acs$day_numb)
daily_acs$ID <- as.numeric(daily_acs$ID)

# transform to df
daily_acs_df <- as.data.frame(daily_acs)
daily_arrest_count_buffer_df <- as.data.frame(daily_arrest_count_buffer)

# merge
daily_arrest_count_buffer_covar <- merge(daily_acs_df, daily_arrest_count_buffer_df, by=c("ID","day_numb", "Name", "Address","City","Zip","Full_Add", "State", "geo_method"))
daily_arrest_count_buffer_covar <- daily_arrest_count_buffer_covar %>% arrange(ID, day_numb)

# drop geometry 
daily_arrest_count_buffer_covar <-  subset(daily_arrest_count_buffer_covar, select = -c(geometry.y,geometry.x))

######################################
#### RUN SYNTHETIC CONTROLS  ##########
######################################

# note: intervention date is Nov 30, 2021: day_numb = 2891; post-period is day_numb 2892 - 3560

# prep data for SCM: create separate df for OPEH and OPWH

# opeh
opeh_ssp_df <- subset(daily_arrest_count_buffer_covar, ID != "2")

# opwh
opwh_ssp_df <- subset(daily_arrest_count_buffer_covar, ID != "1")

######################################
##### ONPOINT: EAST HARLEM  ##########
######################################

# data prep step for synth control
opeh_ssp_dataprep.out<- dataprep(foo = opeh_ssp_df, 
                                           time.predictors.prior = c(1:2891),
                                           predictors=c("median_age",
                                                        "median_household_income_in_past_12_months",
                                                        "pct_black",
                                                        "pct_latino",
                                                        "pct_white",
                                                        "pct_female",
                                                        "pct_male",
                                                        "pct_hs_ged_attainment",
                                                        "pct_poverty",
                                                        "pct_public_assistance",
                                                        "pct_unemployed"),
                                           predictors.op= "mean",
                                           special.predictors = list(list("tot_day_count", 1:2891, "mean")),
                                           dependent =  "tot_day_count", 
                                           unit.variable =  "ID",
                                           unit.names.variable = "Name",
                                           time.variable =  "day_numb",
                                           treatment.identifier = c(1),
                                           controls.identifier = c(3:27),
                                           time.optimize.ssr = c(1:2891),
                                           time.plot = c(1:3560)) #full period Jan 2014 - Sept. 2023
# transform data
opeh_ssp_synth.out = synth(data.prep.obj = opeh_ssp_dataprep.out, method="BFGS",Sigf.ipop = 5)

opeh_ssp_gaps = opeh_ssp_dataprep.out$Y1plot - (opeh_ssp_dataprep.out$Y0plot %*% opeh_ssp_synth.out$solution.w)

opeh_ssp_synth.tables = synth.tab(dataprep.res = opeh_ssp_dataprep.out, synth.res = opeh_ssp_synth.out)


# check gaps 
mean(opeh_ssp_gaps[1:2891]) #before to day-of intervention implementation
mean(opeh_ssp_gaps[2892:3560]) #after intervention

# get average daily arrests before and after intervention: treated site
# before
mean(opeh_ssp_dataprep.out$Y1plot[1:2891])
# after
mean(opeh_ssp_dataprep.out$Y1plot[2892:3560])

# get average daily arrests before and after intervention: control 
# before
mean((opeh_ssp_dataprep.out$Y0plot %*% opeh_ssp_synth.out$solution.w)[1:2891])
# after
mean((opeh_ssp_dataprep.out$Y0plot %*% opeh_ssp_synth.out$solution.w)[2892:3560]) #after, control

######################################
####### PLACEBO TEST: OPEH  ##########
######################################

# prep for parallel processing
n.cores <- parallel::detectCores()-1
my.cluster <- parallel::makeCluster(n.cores, type="PSOCK")
print(my.cluster)
doParallel::registerDoParallel(cl=my.cluster)

# create new ID var to ensure consecutive indicator for looping
opeh_ssp_df$ID2 <- rep(1:27, each=3560)
# original ID = 1; ID2 = 1

# run placebo tests while storing gap data

opeh_ssp_placebo_gaps <- foreach(i=1:27, .packages="Synth", .combine='cbind') %dopar%
  {
    opeh_ssp_placebo_dataprep.out<-
      dataprep(foo = opwh_ssp_df,
               time.predictors.prior = c(1:2891),
               predictors=c("median_age",
                            "median_household_income_in_past_12_months",
                            "pct_black",
                            "pct_latino",
                            "pct_white",
                            "pct_female",
                            "pct_male",
                            "pct_hs_ged_attainment",
                            "pct_poverty",
                            "pct_public_assistance",
                            "pct_unemployed"),
               predictors.op="mean",
               special.predictors = list(list("tot_day_count", 1:2891, "mean")),
               dependent =  "tot_day_count", 
               unit.variable =  "ID2",
               unit.names.variable = "Name",
               time.variable =  "day_numb",
               treatment.identifier = i,
               controls.identifier = (c(1:27)[-i]),
               time.optimize.ssr = c(1:2891),
               time.plot = c(1:3560))
    # run synth
    opeh_ssp_placebo_synth.out <- synth(
      data.prep.obj = opeh_ssp_placebo_dataprep.out,
      method = "BFGS", 
      Sigf.ipop = 5)
    
    # store gaps
    opeh_ssp_placebo_dataprep.out$Y1plot - (opeh_ssp_placebo_dataprep.out$Y0plot %*% opeh_ssp_placebo_synth.out$solution.w)
  }


# make figure
data_opeh_ssp_placebo <- opeh_ssp_placebo_gaps
rownames(data_opeh_ssp_placebo) <- 1:3560

# Set bounds in gaps data
gap.start_opeh_ssp_placebo <- 1
gap.end_opeh_ssp_placebo <- nrow(data_opeh_ssp_placebo)
days <- 1:3560
gap.end.pre_opeh_ssp_placebo  <- which(rownames(data_opeh_ssp_placebo)=="2891")

# Plot
plot(days,data_opeh_ssp_placebo[gap.start_opeh_ssp_placebo:gap.end_opeh_ssp_placebo,which(colnames(data_opeh_ssp_placebo)=="1")],
     main="Placebo Runs: OPEH Daily Arrests",
     ylim=c(-15,20),xlab="Days",
     xlim=c(1,3560),ylab="Gap in drug related arrests",
     type="l",lwd=2,col="black",
     xaxs="i",yaxs="i")

# Add lines for control states
for (i in 1:ncol(data_opeh_ssp_placebo)) { lines(days,data_opeh_ssp_placebo[gap.start_opeh_ssp_placebo:gap.end_opeh_ssp_placebo,i],col="gray") }

# Add treatment Line
lines(days,data_opeh_ssp_placebo[gap.start_opeh_ssp_placebo:gap.end_opeh_ssp_placebo,which(colnames(data_opeh_ssp_placebo)=="1")],lwd=2,col="blue")

# Add grid
abline(v=2892,lty="dotted",lwd=2)
abline(h=0,lty="dashed",lwd=2)

#  R/MSPE Pre-Treatment
all_opeh_ssp_mse_pre <- apply(data_opeh_ssp_placebo[ gap.start_opeh_ssp_placebo:gap.end_opeh_ssp_placebo,]^2,2,mean)
all_opeh_ssp_rmse_pre <- sqrt(all_opeh_ssp_mse_pre)
opeh.mse_pre <- as.numeric(all_opeh_ssp_mse_pre[1])
opeh.rmse_pre <- sqrt(opeh.mse_pre)
print(opeh.rmse_pre)

# R/MSPE post-treatment
all_opeh_ssp_mse_post <- apply(data_opeh_ssp_placebo[(gap.end_opeh_ssp_placebo+1):gap.end_opeh_ssp_placebo,]^2,2,mean)
all_opeh_ssp_rmse_post <- sqrt(all_opeh_ssp_mse_post)

opeh.mse_post <- as.numeric(all_opeh_ssp_mse_post[1])
opeh.rmse_post <- sqrt(opeh.mse_post)
print(opeh.rmse_post)

# Pre/Post
all_opeh_ssp_post_pre_rmspe <- all_opeh_ssp_rmse_post/all_opeh_ssp_rmse_pre
sort(all_opeh_ssp_post_pre_rmspe)

######################################
##### ONPOINT: WASH. HEIGHTS #########
######################################

# data prep step for synth control
opwh_ssp_dataprep.out<- dataprep(foo = opwh_ssp_df, 
                                           time.predictors.prior = c(1:2891),
                                           predictors=c("median_age",
                                                        "median_household_income_in_past_12_months",
                                                        "pct_black",
                                                        "pct_latino",
                                                        "pct_white",
                                                        "pct_female",
                                                        "pct_male",
                                                        "pct_hs_ged_attainment",
                                                        "pct_poverty",
                                                        "pct_public_assistance",
                                                        "pct_unemployed"),
                                           predictors.op= "mean",
                                           special.predictors = list(list("tot_day_count", 1:2891, "mean")),
                                           dependent =  "tot_day_count", 
                                           unit.variable =  "ID",
                                           unit.names.variable = "Name",
                                           time.variable =  "day_numb",
                                           treatment.identifier = c(2),
                                           controls.identifier = c(3:27),
                                           time.optimize.ssr = c(1:2891),
                                           time.plot = c(1:3560)) #full period Jan 2014 - Sept. 2023
#transform data
opwh_ssp_synth.out = synth(data.prep.obj = opwh_ssp_dataprep.out, method="BFGS",Sigf.ipop = 5)

opwh_ssp_gaps = opwh_ssp_dataprep.out$Y1plot - (opwh_ssp_dataprep.out$Y0plot %*% opwh_ssp_synth.out$solution.w)

opwh_ssp_synth.tables = synth.tab(dataprep.res = opwh_ssp_dataprep.out, synth.res = opwh_ssp_synth.out)

# check gaps 
mean(opwh_ssp_gaps[1:2891]) #before treat
mean(opwh_ssp_gaps[2892:3560]) #after treat

# get average daily arrests before and after intervention: treated site
# before
mean(opwh_ssp_dataprep.out$Y1plot[1:2891])
# after
mean(opwh_ssp_dataprep.out$Y1plot[2892:3560])

# get average daily arrests before and after intervention: control 
# before
mean((opwh_ssp_dataprep.out$Y0plot %*% opwh_ssp_synth.out$solution.w)[1:2891])
# after
mean((opwh_ssp_dataprep.out$Y0plot %*% opwh_ssp_synth.out$solution.w)[2892:3560]) #after, control

######################################
####### PLACEBO TEST: OPWH  ##########
######################################

# create new ID var to ensure consecutive indicator for looping
opwh_ssp_df$ID2 <- rep(1:27, each=3560)
# original ID = 1; ID2 = 1

#run placebo tests while storing gap data

opwh_ssp_placebo_gaps <- foreach(i=1:27, .packages="Synth", .combine='cbind') %dopar%
  {
    opwh_ssp_placebo_dataprep.out<-
      dataprep(foo = opwh_ssp_df,
               time.predictors.prior = c(1:2891),
               predictors=c("median_age",
                            "median_household_income_in_past_12_months",
                            "pct_black",
                            "pct_latino",
                            "pct_white",
                            "pct_female",
                            "pct_male",
                            "pct_hs_ged_attainment",
                            "pct_poverty",
                            "pct_public_assistance",
                            "pct_unemployed"),
               predictors.op="mean",
               special.predictors = list(list("tot_day_count", 1:2891, "mean")),
               dependent =  "tot_day_count", 
               unit.variable =  "ID2",
               unit.names.variable = "Name",
               time.variable =  "day_numb",
               treatment.identifier = i,
               controls.identifier = (c(1:27)[-i]),
               time.optimize.ssr = c(1:2891),
               time.plot = c(1:3560))
    # run synth
    opwh_ssp_placebo_synth.out <- synth(
      data.prep.obj = opwh_ssp_placebo_dataprep.out,
      method = "BFGS", 
      Sigf.ipop = 5)
    
    # store gaps
    opwh_ssp_placebo_dataprep.out$Y1plot - (opwh_ssp_placebo_dataprep.out$Y0plot %*% opwh_ssp_placebo_synth.out$solution.w)
  }

# make figure
data_opwh_ssp_placebo <- opwh_ssp_placebo_gaps
rownames(data_opwh_ssp_placebo) <- 1:3560

# Set bounds in gaps data
gap.start_opwh_ssp_placebo <- 1
gap.end_opwh_ssp_placebo <- nrow(data_opwh_ssp_placebo)
gap.end.pre_opwh_ssp_placebo  <- which(rownames(data_opwh_ssp_placebo)=="2891")

# Plot
plot(days,data_opwh_ssp_placebo[gap.start_opwh_ssp_placebo:gap.end_opwh_ssp_placebo,which(colnames(data_opwh_ssp_placebo)=="1")],
     main="Placebo Runs: OPWH Daily Arrests",
     ylim=c(-15,20),xlab="Days",
     xlim=c(1,3560),ylab="Gap in drug related arrests",
     type="l",lwd=2,col="black",
     xaxs="i",yaxs="i")

# Add lines for control states
for (i in 1:ncol(data_opwh_ssp_placebo)) { lines(days,data_opwh_ssp_placebo[gap.start_opwh_ssp_placebo:gap.end_opwh_ssp_placebo,i],col="gray") }

# Add treatment Line
lines(days,data_opwh_ssp_placebo[gap.start_opwh_ssp_placebo:gap.end_opwh_ssp_placebo,which(colnames(data_opwh_ssp_placebo)=="1")],lwd=2,col="blue")

# Add grid
abline(v=2892,lty="dotted",lwd=2)
abline(h=0,lty="dashed",lwd=2)

#  R/MSPE Pre-Treatment
all_opwh_ssp_mse_pre <- apply(data_opwh_ssp_placebo[ gap.start_opwh_ssp_placebo:gap.end_opwh_ssp_placebo,]^2,2,mean)
all_opwh_ssp_rmse_pre <- sqrt(all_opwh_ssp_mse_pre)
opwh.mse_pre <- as.numeric(all_opwh_ssp_mse_pre[1])
opwh.rmse_pre <- sqrt(opwh.mse_pre)
print(opwh.rmse_pre)

# R/MSPE post-treatment
all_opwh_ssp_mse_post <- apply(data_opwh_ssp_placebo[(gap.end_opwh_ssp_placebo+1):gap.end_opwh_ssp_placebo,]^2,2,mean)
all_opwh_ssp_rmse_post <- sqrt(all_opwh_ssp_mse_post)

opwh.mse_post <- as.numeric(all_opwh_ssp_mse_post[1])
opwh.rmse_post <- sqrt(opwh.mse_post)
print(opwh.rmse_post)

# Pre/Post
all_opwh_ssp_post_pre_rmspe <- all_opwh_ssp_rmse_post/all_opwh_ssp_rmse_pre
sort(all_opwh_ssp_post_pre_rmspe)
