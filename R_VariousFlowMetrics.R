library(tidyverse)
library(lubridate)
library(AICcmodavg)
library(ggpubr)
library(waterYearType)

# DOWNLOAD wy 1997-2020 FROM DAYFLOW  ---------------------------------------------

wy_1997_2020 <- read.csv('https://data.cnra.ca.gov/dataset/06ee2016-b138-47d7-9e85-f46fae674536/resource/21c377fe-53b8-4bd6-9e1f-2025221be095/download/dayflow-results-1997-2020.csv',
                         stringsAsFactors = FALSE)
wy_1955_1996 <- read.csv("dayflowCalculations_pre97.csv")

dayflow <- bind_rows(wy_1955_1996,wy_1997_2020)

remove(wy_1955_1996,wy_1997_2020)

dayflow$Date<-as.Date(dayflow$Date,"%m/%d/%Y")
dayflow$Month<-month(dayflow$Date)

#Export out in case we need it later
monthly_dayflow <- dayflow %>% 
  group_by(Year, Month) %>% 
  summarise(X2 = mean(X2, na.rm = TRUE),
            Outflow=mean(OUT, na.rm = TRUE),
            Inflow=mean(TOT, na.rm = TRUE))

write.csv(monthly_dayflow, file = 'dayflow_monthly.csv', row.names = FALSE)


#Show monthly inflow vs monthly outflow linear regression results
summary(lm(monthly_dayflow$Inflow~monthly_dayflow$Outflow))
#Multiple R-squared:  0.9886,	Adjusted R-squared:  0.9886 


#Species to analyze: Longfin Smelt (FMWT), Splittail (FMWT), Starry Flounder (Bay Study), Bay Shrimp (Bay Study)

#Create function
longfin_flow<-function(variable){
  mean(variable[Month %in% 1:6])
}


#Spring flow summarized as in Kimmerer 2002 table 4
spring_flow <- dayflow %>% 
  group_by(Year) %>%  
  summarise(X2_janjun = mean(X2[Month %in% 1:6]), # LFS
            X2_marjun = mean(X2[Month %in% 3:6]), # Starry Flounder
            X2_febmay = mean(X2[Month %in% 2:5]), # Splittail
            X2_marmay = mean(X2[Month %in% 3:5]), # Bay Shrimp
            X2_junjul = mean(X2[Month %in% 6:7]), # Striped Bass, see Tamburello et al.
            Outflow_janjun = mean(OUT[Month %in% 1:6]), # LFS
            Outflow_marjun = mean(OUT[Month %in% 3:6]), # Starry Flounder
            Outflow_febmay = mean(OUT[Month %in% 2:5]), # Splittail
            Outflow_marmay = mean(OUT[Month %in% 3:5]), # Bay Shrimp
            Outflow_junjul = mean(OUT[Month %in% 6:7]), # Striped Bass
            Inflow_janjun = mean(TOT[Month %in% 1:6]), # LFS
            Inflow_marjun = mean(TOT[Month %in% 3:6]), # Starry Flounder
            Inflow_febmay = mean(TOT[Month %in% 2:5]), # Splittail
            Inflow_marmay = mean(TOT[Month %in% 3:5]), # Bay Shrimp
            Inflow_junjul = mean(TOT[Month %in% 6:7]), # Striped Bass
  )

#Starry flounder is supposed to use lag (previous year metric)
spring_flow$X2_marjun_prev_year<-lag(spring_flow$X2_marjun)
spring_flow$Outflow_marjun_prev_year<-lag(spring_flow$Outflow_marjun)
spring_flow$Inflow_marjun_prev_year<-lag(spring_flow$Inflow_marjun)

#Grab WSIHIST data for unimpaired runoff
WSIHIST<- water_year_indices %>% rename(Year=WY) %>% group_by(Year) %>% summarise(WYsum=sum(WYsum))

full_flow_data<-left_join(spring_flow,WSIHIST) %>% rename(Unimpaired_runoff=WYsum)

##############################################################################
#Load biological data
master_data <- read.csv('MetWater_Database_111017_updated_032521_R.csv')

#Merge with predictor variables
master_data<-left_join(master_data,full_flow_data)

#create another categorical dummy variable for pre-clam, post-clam, and post-POD
master_data$step_three<-as.factor(ifelse(master_data$Year<1988,"PreClam",ifelse(master_data$Year<2002,"PostClam","PostPOD")))
master_data$step_three<-ordered(master_data$step_three, levels = c("PreClam", "PostClam", "PostPOD"))

#Create the proper abundance values used in Kimmerer papers

#Species:
#Longfin Smelt = fmwt_lfs
#Splittail = fmwt_ssa
#Starry Flounder = bsot_sf
#Bay Shrimp = bsot_bsh
#Striped Bass = tns_log_sb
master_data<- master_data %>% mutate(log_lfs = log10(fmwt_lfs),log_spt = log10(fmwt_ssa + 1),log_stf = log10(bsot_sf + 10),log_bsh = log10(bsot_bsh))

#Add shorthand for year
master_data$year_short<-sprintf('%02d',master_data$Year %% 100)



################################################################################
#Redo Kimmerer et al. 2009 models - with updated data
################################################################################

linearmodel_lfs_updated_x2<- lm(log_lfs ~ X2_janjun + step_1987, data = master_data)
linearmodel_stb_updated_x2<- lm(tns_log_sb ~ X2_junjul, data = master_data)
linearmodel_splt_updated_x2<- lm(log_spt ~ X2_febmay , data = master_data)

linearmodel_lfs_updated_outflow<- lm(log_lfs ~ Outflow_janjun + step_1987, data = master_data)
linearmodel_stb_updated_outflow<- lm(tns_log_sb ~ Outflow_junjul, data = master_data)
linearmodel_splt_updated_outflow<- lm(log_spt ~ Outflow_febmay , data = master_data)

linearmodel_lfs_updated_inflow<- lm(log_lfs ~ Inflow_janjun + step_1987, data = master_data)
linearmodel_stb_updated_inflow<- lm(tns_log_sb ~ Inflow_junjul, data = master_data)
linearmodel_splt_updated_inflow<- lm(log_spt ~ Inflow_febmay , data = master_data)

linearmodel_lfs_updated_Unimpaired_runoff<- lm(log_lfs ~ Unimpaired_runoff + step_1987, data = master_data)
linearmodel_stb_updated_Unimpaired_runoff<- lm(tns_log_sb ~Unimpaired_runoff, data = master_data)
linearmodel_splt_updated_Unimpaired_runoff<- lm(log_spt ~ Unimpaired_runoff , data = master_data)

sum_table<-data.frame(species_relationship = c('Longfin Smelt','Striped Bass','Splittail'),
                      X2= c(summary(linearmodel_lfs_updated_x2)$r.squared,summary(linearmodel_stb_updated_x2)$r.squared,summary(linearmodel_splt_updated_x2)$r.squared),
                      Outflow = c(summary(linearmodel_lfs_updated_outflow)$r.squared,summary(linearmodel_stb_updated_outflow)$r.squared,summary(linearmodel_splt_updated_outflow)$r.squared),
                      Inflow = c(summary(linearmodel_lfs_updated_inflow)$r.squared,summary(linearmodel_stb_updated_inflow)$r.squared,summary(linearmodel_splt_updated_inflow)$r.squared),
                      Unimpaired_runoff = c(summary(linearmodel_lfs_updated_Unimpaired_runoff)$r.squared,summary(linearmodel_stb_updated_Unimpaired_runoff)$r.squared,summary(linearmodel_splt_updated_Unimpaired_runoff)$r.squared)
)

write.csv(sum_table,file="Summary_R-squared_Table_for_flow_metrics.csv",row.names = F)
