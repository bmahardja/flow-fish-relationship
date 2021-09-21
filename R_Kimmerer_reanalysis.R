library(tidyverse)
library(lubridate)
library(AICcmodavg)
library(ggpubr)

setwd("D:/2021-03-25 - Kimmerer reanalysis memo")


###########################################################################################################################
#Set up data
###########################################################################################################################

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
  summarise(X2 = mean(X2, na.rm = TRUE))

write.csv(monthly_dayflow, file = 'dayflow_monthly.csv', row.names = FALSE)

#Species to analyze: Longfin Smelt (FMWT), Splittail (FMWT), Starry Flounder (Bay Study), Bay Shrimp (Bay Study)

#Spring flow summarized as in Kimmerer 2002 table 4
spring_flow <- dayflow %>% 
  group_by(Year) %>%  
  summarise(X2_janjun = mean(X2[Month %in% 1:6]), # LFS
            X2_marjun = mean(X2[Month %in% 3:6]), # Starry Flounder
            X2_febmay = mean(X2[Month %in% 2:5]), # Splittail
            X2_marmay = mean(X2[Month %in% 3:5]) # Bay Shrimp
  )

#Starry flounder is supposed to use lag (previous year metric)
spring_flow$X2_marjun_prev_year<-lag(spring_flow$X2_marjun)

#Load biological data
master_data <- read.csv('MetWater_Database_111017_updated_032521_R.csv')

#Merge with predictor variables
master_data<-left_join(master_data,spring_flow)

#create another categorical dummy variable for pre-clam, post-clam, and post-POD
master_data$step_three<-as.factor(ifelse(master_data$Year<1988,"PreClam",ifelse(master_data$Year<2002,"PostClam","PostPOD")))
master_data$step_three<-ordered(master_data$step_three, levels = c("PreClam", "PostClam", "PostPOD"))

#Create the proper abundance values used in Kimmerer papers

#Species:
#Longfin Smelt = fmwt_lfs
#Splittail = fmwt_ssa
#Starry Flounder = bsot_sf
#Bay Shrimp = bsot_bsh
master_data<- master_data %>% mutate(log_lfs = log10(fmwt_lfs),log_spt = log10(fmwt_ssa + 1),log_stf = log10(bsot_sf + 10),log_bsh = log10(bsot_bsh))

#Add shorthand for year
master_data$year_short<-sprintf('%02d',master_data$Year %% 100)

###########################################################################################################################
#Re-create Kimmerer et al. 2009 models
###########################################################################################################################

kimmerer2009_lfs<- lm(log_lfs ~ X2_janjun + step_1987, data = (master_data %>% filter(Year<=2007)))
kimmerer2009_splt<- lm(log_spt ~ X2_febmay , data = (master_data %>% filter(Year<=2007)))
kimmerer2009_stf<- lm(log_stf ~ X2_marjun_prev_year + step_1987, data = (master_data %>% filter(Year<=2007)))
#Note that paper states that bay shrimp only until 2006
kimmerer2009_bsh<- lm(log_bsh ~ X2_marmay, data = (master_data %>% filter(Year<=2006)))

#Create function to create data frame
create_data<-function(df,model){
  new_data<-df %>% mutate(Prediction=predict(model)) %>% 
    left_join(master_data) %>% 
    select(c(1:4,"Year","year_short"))
  return(new_data)
}

test<-create_data(df=kimmerer2009_bsh$model,model=kimmerer2009_bsh)


#Add different color combination
color_set<-c("#4363d8","#f58231","#800000")

#Create function to create ggplot for each model
generic_plot<-function(df,fig_title,x_axis,y_axis,rsquared){
  #For step change at 1987
  if("step_1987" %in% colnames(df) ){
    newplot<-ggplot(data=df)+
      theme_bw()+
      geom_point(size=7,aes(x=df[,2],y=df[,1], colour=as.factor(df[,3])),alpha=0.6)+
      geom_line(aes(x=df[,2],y = df[,4], colour=as.factor(df[,3])),size=1.5,linetype="dashed")+ 
      geom_text(aes(x=df[,2],y=df[,1],label=year_short),hjust=0,vjust=0,size=6,alpha=0.5)+
      theme(plot.title = element_text(size=28),axis.text.x = element_text(size=21, color="black"),
            axis.text.y = element_text(size=20, color="black"),axis.title.x = element_text(size = 27, angle = 00),
            axis.title.y = element_text(size = 27, angle = 90),legend.title = element_blank(),
            legend.position = c(0.85, 0.85),legend.text = element_text(face="bold",size=20),
            legend.background = element_rect(size=1, linetype="solid",colour="black"))+
      #scale_y_continuous(limits = c(0, 5))+
      scale_color_manual(values=color_set, labels = c("pre-clam","post-clam"))+
      labs(title=paste(fig_title,", R^2 = ",rsquared,sep=""),x=x_axis, y=y_axis)
    return(newplot) 
  }
  #For 3 different eras
  if("step_three" %in% colnames(df)){
    newplot<-ggplot(data=df)+
      theme_bw()+
      geom_point(size=7,aes(x=df[,2],y=df[,1], colour=df[,3]),alpha=0.6)+
      geom_line(aes(x=df[,2],y = df[,4], colour=df[,3]),size=1.5,linetype="dashed")+ 
      geom_text(aes(x=df[,2],y=df[,1],label=year_short),hjust=0,vjust=0,size=6,alpha=0.5)+
      theme(plot.title = element_text(size=28),axis.text.x = element_text(size=21, color="black"),
            axis.text.y = element_text(size=20, color="black"),axis.title.x = element_text(size = 27, angle = 00),
            axis.title.y = element_text(size = 27, angle = 90),legend.title = element_blank(),
            legend.position = c(0.85, 0.85),legend.text = element_text(face="bold",size=20),
            legend.background = element_rect(size=1, linetype="solid",colour="black"))+
      #scale_y_continuous(limits = c(0, 5))+
      scale_color_manual(values=color_set, labels = c("pre-clam","post-clam","post-POD"))+
      labs(title=paste(fig_title,", R^2 = ",rsquared,sep=""),x=x_axis, y=y_axis)
    return(newplot)
  }
  #For no dummy variables
  if(!("step_1987" %in% colnames(data)&"step_three" %in% colnames(df))){
    newplot<-ggplot(data=df %>% mutate(step_1987=as.factor(ifelse(Year<1988,0,1))))+
      theme_bw()+
      geom_point(size=7,aes(x=df[,2],y=df[,1], colour=step_1987),alpha=0.6)+
      geom_line(aes(x=df[,2],y = df[,3]),size=1.5,linetype="dashed")+ 
      geom_text(aes(x=df[,2],y=df[,1],label=year_short),hjust=0,vjust=0,size=6,alpha=0.5)+
      theme(plot.title = element_text(size=28),axis.text.x = element_text(size=21, color="black"),
            axis.text.y = element_text(size=20, color="black"),axis.title.x = element_text(size = 27, angle = 00),
            axis.title.y = element_text(size = 27, angle = 90),legend.title = element_blank(),
            legend.position = c(0.85, 0.85),legend.text = element_text(face="bold",size=20),
            legend.background = element_rect(size=1, linetype="solid",colour="black"))+
      #scale_y_continuous(limits = c(0, 5))+
      scale_color_manual(values=color_set, labels = c("pre-clam","post-clam"))+
      #annotate('text', x = 0, y = 0, label = expression(paste(R^{2},'=',rsquared),parse = TRUE,size=20)+
      labs(title=paste(fig_title,", R^2 = ",rsquared,sep=""),x=x_axis, y=y_axis)
    return(newplot)
  }
}

#Plots for Kimmerer et al. 2009 models
plot_kimmerer2009_lfs<-generic_plot(df=create_data(df=kimmerer2009_lfs$model,model=kimmerer2009_lfs),rsquared=round(summary(kimmerer2009_bsh)$adj.r.squared,digits=2),fig_title = "Longfin Smelt",x_axis="X2",y_axis = "Log Abundance")
plot_kimmerer2009_bsh<-generic_plot(df=create_data(df=kimmerer2009_bsh$model,model=kimmerer2009_bsh),rsquared=round(summary(kimmerer2009_bsh)$adj.r.squared,digits=2),fig_title = "Bay Shrimp",x_axis="X2",y_axis = "Log Abundance")
plot_kimmerer2009_stf<-generic_plot(df=create_data(df=kimmerer2009_stf$model,model=kimmerer2009_stf),rsquared=round(summary(kimmerer2009_bsh)$adj.r.squared,digits=2),fig_title = "Starry Flounder",x_axis="X2",y_axis = "Log Abundance + 10")
plot_kimmerer2009_splt<-generic_plot(df=create_data(df=kimmerer2009_splt$model,model=kimmerer2009_splt),rsquared=round(summary(kimmerer2009_bsh)$adj.r.squared,digits=2),fig_title = "Splittail",x_axis="X2",y_axis = "Log Abundance + 1")

#plot_kimmerer2009_bsh
#plot_kimmerer2009_lfs
#plot_kimmerer2009_stf
#plot_kimmerer2009_splt

summary(kimmerer2009_splt)$adj.r.squared

################################################################################
#Redo Kimmerer et al. 2009 models - with updated data
################################################################################

kimmerer2009_lfs_updated<- lm(log_lfs ~ X2_janjun + step_1987, data = master_data)
kimmerer2009_splt_updated<- lm(log_spt ~ X2_febmay , data = master_data)
kimmerer2009_stf_updated<- lm(log_stf ~ X2_marjun_prev_year + step_1987, data = master_data)
kimmerer2009_bsh_updated<- lm(log_bsh ~ X2_marmay, data = master_data)
summary(kimmerer2009_lfs_updated)

plot_kimmerer2009_lfs_updated<-generic_plot(df=create_data(df=kimmerer2009_lfs_updated$model,model=kimmerer2009_lfs_updated),rsquared=round(summary(kimmerer2009_lfs_updated)$adj.r.squared,digits=2),fig_title = "Longfin Smelt",x_axis="X2",y_axis = "Log Abundance")
plot_kimmerer2009_bsh_updated<-generic_plot(df=create_data(df=kimmerer2009_bsh_updated$model,model=kimmerer2009_bsh_updated),rsquared=round(summary(kimmerer2009_bsh_updated)$adj.r.squared,digits=2),fig_title = "Bay Shrimp",x_axis="X2",y_axis = "Log Abundance")
plot_kimmerer2009_stf_updated<-generic_plot(df=create_data(df=kimmerer2009_stf_updated$model,model=kimmerer2009_stf_updated),rsquared=round(summary(kimmerer2009_stf_updated)$adj.r.squared,digits=2),fig_title = "Starry Flounder",x_axis="X2",y_axis = "Log Abundance + 10")
plot_kimmerer2009_splt_updated<-generic_plot(df=create_data(df=kimmerer2009_splt_updated$model,model=kimmerer2009_splt_updated),rsquared=round(summary(kimmerer2009_splt_updated)$adj.r.squared,digits=2),fig_title = "Splittail",x_axis="X2",y_axis = "Log Abundance + 1")

#plot_kimmerer2009_bsh_updated
#plot_kimmerer2009_lfs_updated
#plot_kimmerer2009_stf_updated
#plot_kimmerer2009_splt_updated

plot_kimmerer2009_summary<-ggarrange(ggarrange(plot_kimmerer2009_lfs, plot_kimmerer2009_splt, plot_kimmerer2009_stf, plot_kimmerer2009_bsh, ncol=1, nrow=4), 
                                     ggarrange(plot_kimmerer2009_lfs_updated, plot_kimmerer2009_splt_updated, plot_kimmerer2009_stf_updated, plot_kimmerer2009_bsh_updated, ncol=1, nrow=4), ncol=2, nrow=1)
  
plot_kimmerer2009_summary<-annotate_figure(plot_kimmerer2009_summary,top=text_grob("Kimmerer et al. 2009 original models                                         Models updated with more recent data",size=28,x=0.5))

#Print out original and updated figures side-by-side
png(filename="Kimmerer et al. 2009 updated.png", units="in",type="cairo", bg="white", height=36, 
      width=24, res=300, pointsize=20)
plot_kimmerer2009_summary
dev.off()


##################################################
##Summary Error Bars figure
##################################################

#Create function to pull X2 slope coefficient and 95% conf int
create_confint<-function(model,model.type,species.name){
  new_data<-confint(model, level=0.95) %>% as_tibble() %>% slice(2) %>%
    mutate(coefficient=model$coefficient[2],model_type=model.type,species=species.name) %>% 
    rename(L95='2.5 %',U95='97.5 %') 
  return(new_data)
}


#test<-create_confint(model=kimmerer2009_lfs_updated,model.type="Updated",species.name="Longfin Smelt")

data_model_sum<-bind_rows(create_confint(model=kimmerer2009_lfs,model.type="Original",species.name="Longfin Smelt"),
                          create_confint(model=kimmerer2009_splt,model.type="Original",species.name="Splittail"),
                          create_confint(model=kimmerer2009_stf,model.type="Original",species.name="Starry Flounder"),
                          create_confint(model=kimmerer2009_bsh,model.type="Original",species.name="Bay Shrimp"),
                          create_confint(model=kimmerer2009_lfs_updated,model.type="Updated",species.name="Longfin Smelt"),
                          create_confint(model=kimmerer2009_splt_updated,model.type="Updated",species.name="Splittail"),
                          create_confint(model=kimmerer2009_stf_updated,model.type="Updated",species.name="Starry Flounder"),
                          create_confint(model=kimmerer2009_bsh_updated,model.type="Updated",species.name="Bay Shrimp")
)


data_model_sum$model_type<-factor(data_model_sum$model_type)
data_model_sum$species<-factor(data_model_sum$species)
data_model_sum

#Create figure to summarize changes to slope coefficient
plot_error_bar <- ggplot(data_model_sum, aes(x=model_type, y=coefficient)) + 
  facet_grid(cols = vars(species))+
  geom_line(aes(group = species), lty = 2)+ geom_point(size=4)+
  geom_errorbar(aes(ymin=L95, ymax=U95),width=0.5)+
  theme(plot.title = element_text(size=28),axis.text.x = element_text(size=14, color="black"),
        axis.text.y = element_text(size=18, color="black"),axis.title.x = element_text(size = 22, angle = 00),
        axis.title.y = element_text(size = 22, angle = 90),strip.text.x = element_text(size = 18))+
  labs(x="Dataset", y="X2 coefficient estimate")

plot_error_bar


#Print out figure
png(filename="Error bar plot.png", units="in",type="cairo", bg="white", height=8, 
    width=16, res=300, pointsize=20)
plot_error_bar
dev.off()





##################################################
##Model selection
##################################################


#Longfin Smelt
Cand.set.LFS <- list( )
Cand.set.LFS[[1]] <- lm(log_lfs ~ X2_janjun, data = master_data)
Cand.set.LFS[[2]] <- lm(log_lfs ~ X2_janjun + step_1987, data = master_data)
Cand.set.LFS[[3]] <- lm(log_lfs ~ X2_janjun + step_1987 + X2_janjun:step_1987, data = master_data)
Cand.set.LFS[[4]] <- lm(log_lfs ~ X2_janjun + step_three, data = master_data)
Cand.set.LFS[[5]] <- lm(log_lfs ~ X2_janjun + step_three + X2_janjun:step_three, data = master_data)

##create a vector of names to trace back models in set
Modnames.LFS <- paste("mod","LFS", 1:length(Cand.set.LFS), sep = "_")
##generate AICc table
aictab(cand.set = Cand.set.LFS, modnames = Modnames.LFS, sort = TRUE)
table_LFS<-aictab(cand.set = Cand.set.LFS, modnames = Modnames.LFS, sort = TRUE)

#Model 4 was best by 0.86 weight
model_best_lfs <- lm(log_lfs ~ X2_janjun + step_three, data = master_data)
#Create figure
plot_model_best_lfs<-generic_plot(df=create_data(df=model_best_lfs$model,model=model_best_lfs),rsquared=round(summary(model_best_lfs)$adj.r.squared,digits=2),fig_title = "Longfin Smelt",x_axis="X2",y_axis = "Log Abundance")
plot_model_best_lfs



#Splittail
Cand.set.SPLT <- list( )
Cand.set.SPLT[[1]] <- lm(log_spt ~ X2_febmay, data = master_data)
Cand.set.SPLT[[2]] <- lm(log_spt ~ X2_febmay + step_1987, data = master_data)
Cand.set.SPLT[[3]] <- lm(log_spt ~ X2_febmay + step_1987 + X2_febmay:step_1987, data = master_data)
Cand.set.SPLT[[4]] <- lm(log_spt ~ X2_febmay + step_three, data = master_data)
Cand.set.SPLT[[5]] <- lm(log_spt ~ X2_febmay + step_three + X2_febmay:step_three, data = master_data)
##create a vector of names to trace back models in set
Modnames.SPLT <- paste("mod","SPLT", 1:length(Cand.set.SPLT), sep = "_")
##generate AICc table
aictab(cand.set = Cand.set.SPLT, modnames = Modnames.SPLT, sort = TRUE)
table_SPLT<-aictab(cand.set = Cand.set.SPLT, modnames = Modnames.SPLT, sort = TRUE)

#Model 4 was best by 0.81 weight
model_best_splt <- lm(log_spt ~ X2_febmay + step_three, data = master_data)
summary(model_best_splt)
#Create figure
plot_model_best_splt<-generic_plot(df=create_data(df=model_best_splt$model,model=model_best_splt),rsquared=round(summary(model_best_splt)$adj.r.squared,digits=2),fig_title = "Splittail",x_axis="X2",y_axis = "Log Abundance + 1")
plot_model_best_splt



#Starry Flounder
Cand.set.STF <- list( )
Cand.set.STF[[1]] <- lm(log_stf ~ X2_marjun_prev_year, data = master_data)
Cand.set.STF[[2]] <- lm(log_stf ~ X2_marjun_prev_year + step_1987, data = master_data)
Cand.set.STF[[3]] <- lm(log_stf ~ X2_marjun_prev_year + step_1987 + X2_marjun_prev_year:step_1987, data = master_data)
Cand.set.STF[[4]] <- lm(log_stf ~ X2_marjun_prev_year + step_three, data = master_data)
Cand.set.STF[[5]] <- lm(log_stf ~ X2_marjun_prev_year + step_three + X2_marjun_prev_year:step_three, data = master_data)
##create a vector of names to trace back models in set
Modnames.STF <- paste("mod","STF", 1:length(Cand.set.STF), sep = "_")
##generate AICc table
aictab(cand.set = Cand.set.STF, modnames = Modnames.STF, sort = TRUE)
table_STF<-aictab(cand.set = Cand.set.STF, modnames = Modnames.STF, sort = TRUE)

#Model 2 was best by 0.64 weight
model_best_stf <- lm(log_stf ~ X2_marjun_prev_year + step_1987, data = master_data)
#Create figure
plot_model_best_stf<-generic_plot(df=create_data(df=model_best_stf$model,model=model_best_stf),rsquared=round(summary(model_best_stf)$adj.r.squared,digits=2),fig_title = "Starry Flounder",x_axis="X2",y_axis = "Log Abundance + 10")
plot_model_best_stf


#Bay Shrimp
Cand.set.BSH <- list( )
Cand.set.BSH[[1]] <- lm(log_bsh ~ X2_marmay, data = master_data)
Cand.set.BSH[[2]] <- lm(log_bsh ~ X2_marmay + step_1987, data = master_data)
Cand.set.BSH[[3]] <- lm(log_bsh ~ X2_marmay + step_1987 + X2_marmay:step_1987, data = master_data)
Cand.set.BSH[[4]] <- lm(log_bsh ~ X2_marmay + step_three, data = master_data)
Cand.set.BSH[[5]] <- lm(log_bsh ~ X2_marmay + step_three + X2_marmay:step_three, data = master_data)
##create a vector of names to trace back models in set
Modnames.BSH <- paste("mod","BSH", 1:length(Cand.set.BSH), sep = "_")
##generate AICc table
aictab(cand.set = Cand.set.BSH, modnames = Modnames.BSH, sort = TRUE)
table_BSH<-aictab(cand.set = Cand.set.BSH, modnames = Modnames.BSH, sort = TRUE)

#Model 3 was best by 0.68 weight
model_best_bsh <- lm(log_bsh ~ X2_marmay, data = master_data)
#Create figure
plot_model_best_bsh<-generic_plot(df=create_data(df=model_best_bsh$model,model=model_best_bsh),rsquared=round(summary(model_best_bsh)$adj.r.squared,digits=2),fig_title = "Bay Shrimp",x_axis="X2",y_axis = "Log Abundance")
plot_model_best_bsh




#Print out best models
png(filename="Best model by AICc.png", units="in",type="cairo", bg="white", height=36, 
    width=12, res=300, pointsize=20)
ggarrange(plot_model_best_lfs, plot_model_best_splt, plot_model_best_stf, plot_model_best_bsh, ncol=1, nrow=4) 
dev.off()

#Export tables out
table_sum<-bind_rows(table_LFS,table_SPLT,table_STF,table_BSH)
write.csv(table_sum,file="AICc_Summary_Table.csv",row.names = F)



###########################################################################
#Prediction interval for species for 2 points of X2
###########################################################################

#Comparing 75 vs 85 km
#Create dataset

example_x2=c(75,85,75,85)

example_flow <- data.frame(X2_janjun =example_x2, # LFS
                           X2_marjun_prev_year = example_x2, # Starry Flounder
                           X2_febmay = example_x2, # Splittail
                           X2_marmay = example_x2, # Bay Shrimp
                           step_1987 = 1,
                           step_three = c("PostClam","PostClam","PostPOD","PostPOD")
  )


prediction_LFS<-predict(model_best_lfs, example_flow, interval="predict")
prediction_SPLT<-predict(model_best_splt, example_flow, interval="predict")
prediction_STF<-predict(model_best_stf, example_flow, interval="predict")
prediction_BSH<-predict(model_best_bsh, example_flow, interval="predict")

example_flow <- cbind(example_flow,prediction_LFS) %>% rename(LFS_fit=fit,LFS_lwr=lwr,LFS_upr=upr) %>% mutate(LFS_fit=10^(LFS_fit),LFS_lwr=10^(LFS_lwr),LFS_upr=10^(LFS_upr))
example_flow <- cbind(example_flow,prediction_SPLT) %>% rename(SPLT_fit=fit,SPLT_lwr=lwr,SPLT_upr=upr) %>% mutate(SPLT_fit=10^(SPLT_fit)-1,SPLT_lwr=10^(SPLT_lwr)-1,SPLT_upr=10^(SPLT_upr)-1)
example_flow <- cbind(example_flow,prediction_STF) %>% rename(STF_fit=fit,STF_lwr=lwr,STF_upr=upr) %>% mutate(STF_fit=10^(STF_fit)-10,STF_lwr=10^(STF_lwr)-10,STF_upr=10^(STF_upr)-10)
example_flow <- cbind(example_flow,prediction_BSH) %>% rename(BSH_fit=fit,BSH_lwr=lwr,BSH_upr=upr) %>% mutate(BSH_fit=10^(BSH_fit),BSH_lwr=10^(BSH_lwr),BSH_upr=10^(BSH_upr))

#Export table out
write.csv(example_flow,file="Prediction_summary_table.csv",row.names = F)


###########################################################################
#In backtransform scale

mas
example_x2=seq(45,90,0.5)

example_flow <- data.frame(X2_janjun =example_x2, # LFS
                           X2_marjun_prev_year = example_x2, # Starry Flounder
                           X2_febmay = example_x2, # Splittail
                           X2_marmay = example_x2, # Bay Shrimp
                           step_1987 = 1,
                           step_three = "PostPOD"
)



prediction_LFS<-predict(model_best_lfs, example_flow, interval="predict")
prediction_SPLT<-predict(model_best_splt, example_flow, interval="predict")
prediction_STF<-predict(model_best_stf, example_flow, interval="predict")
prediction_BSH<-predict(model_best_bsh, example_flow, interval="predict")

example_flow <- cbind(example_flow,prediction_LFS) %>% rename(LFS_fit=fit,LFS_lwr=lwr,LFS_upr=upr) %>% mutate(LFS_fit=10^(LFS_fit),LFS_lwr=10^(LFS_lwr),LFS_upr=10^(LFS_upr))
example_flow <- cbind(example_flow,prediction_SPLT) %>% rename(SPLT_fit=fit,SPLT_lwr=lwr,SPLT_upr=upr) %>% mutate(SPLT_fit=10^(SPLT_fit)-1,SPLT_lwr=10^(SPLT_lwr)-1,SPLT_upr=10^(SPLT_upr)-1)
example_flow <- cbind(example_flow,prediction_STF) %>% rename(STF_fit=fit,STF_lwr=lwr,STF_upr=upr) %>% mutate(STF_fit=10^(STF_fit)-10,STF_lwr=10^(STF_lwr)-10,STF_upr=10^(STF_upr)-10)
example_flow <- cbind(example_flow,prediction_BSH) %>% rename(BSH_fit=fit,BSH_lwr=lwr,BSH_upr=upr) %>% mutate(BSH_fit=10^(BSH_fit),BSH_lwr=10^(BSH_lwr),BSH_upr=10^(BSH_upr))

data_lfs<-create_data(df=model_best_lfs$model,model=model_best_lfs)
data_lfs$LFS_abundance<-10^(data_lfs$log_lfs)

plot_backtransform_lfs<-ggplot()+
  theme_bw()+
  geom_point(data=data_lfs %>%filter(step_three=="PostPOD"),size=7,aes(x=X2_janjun,y=LFS_abundance),alpha=0.6)+
  geom_line(data=example_flow,aes(x=X2_janjun,y =LFS_fit),size=1.5,linetype="dashed")+ 
  geom_ribbon(data=example_flow,aes(x=X2_janjun, y=LFS_fit, ymax=LFS_upr, ymin=LFS_lwr), alpha=0.2)+
  theme(plot.title = element_text(size=28),axis.text.x = element_text(size=21, color="black"),
        axis.text.y = element_text(size=20, color="black"),axis.title.x = element_text(size = 27, angle = 00),
        axis.title.y = element_text(size = 27, angle = 90),legend.title = element_blank(),
        legend.position = c(0.85, 0.85),legend.text = element_text(face="bold",size=20),
        legend.background = element_rect(size=1, linetype="solid",colour="black"))+
  #scale_y_continuous(limits = c(0, 2000))+
  labs(x="X2", y="Longfin Smelt Abundance Index")


plot_backtransform_lfs


plot_backtransform_lfs_02<-ggplot()+
  theme_bw()+
  geom_point(data=data_lfs %>%filter(step_three=="PostPOD"),size=7,aes(x=X2_janjun,y=LFS_abundance),alpha=0.6)+
  geom_line(data=example_flow,aes(x=X2_janjun,y =LFS_fit),size=1.5,linetype="dashed")+ 
  #geom_ribbon(data=example_flow,aes(x=X2_janjun, y=LFS_fit, ymax=LFS_upr, ymin=LFS_lwr), alpha=0.2)+
  theme(plot.title = element_text(size=28),axis.text.x = element_text(size=21, color="black"),
        axis.text.y = element_text(size=20, color="black"),axis.title.x = element_text(size = 27, angle = 00),
        axis.title.y = element_text(size = 27, angle = 90),legend.title = element_blank(),
        legend.position = c(0.85, 0.85),legend.text = element_text(face="bold",size=20),
        legend.background = element_rect(size=1, linetype="solid",colour="black"))+
  scale_y_continuous(limits = c(0, 2000))+
  geom_text(data=data_lfs %>%filter(step_three=="PostPOD"),aes(x=X2_janjun,y=LFS_abundance,label=year_short),hjust=0,vjust=0,size=6,alpha=0.5)+
  labs(x="X2", y="Longfin Smelt Abundance Index")


plot_backtransform_lfs_02


#Print out figure
png(filename="LFS_backtransformed_plot_01.png", units="in",type="cairo", bg="white", height=12, 
    width=16, res=300, pointsize=20)
plot_backtransform_lfs
dev.off()

png(filename="LFS_backtransformed_plot_02.png", units="in",type="cairo", bg="white", height=12, 
    width=16, res=300, pointsize=20)
plot_backtransform_lfs_02
dev.off()


###########################################################################
#ordered(master_data$step_three, levels = c("PreClam", "PostClam", "PostPOD"))

df<-create_data(df=kimmerer2009_stf$model,model=kimmerer2009_stf)

newplot<-ggplot(data=df)+
  geom_point(size=7,aes(x=df[,2],y=df[,1], colour=as.factor(df[,3])),alpha=0.6)+
  geom_line(aes(x=df[,2],y = df[,4], colour=as.factor(df[,3])),size=1.5,linetype="dashed")+ 
  theme(plot.title = element_text(size=28),axis.text.x = element_text(size=21, color="black"),
        axis.text.y = element_text(size=20, color="black"),axis.title.x = element_text(size = 27, angle = 00),
        axis.title.y = element_text(size = 27, angle = 90))

newplot
df[4]

# SWRCB models fit to most recent data

# Longfin Smelt FMWT ~ January - June X2 + post-2002 Step 
lfs_lm <- lm(log_lfs ~ X2_janjun + step_1987, data = master_data)
summary(lfs_lm)

# Sacramento Splittail FMWT + 1 ~ February - May X2 
spt_lm <- lm(log_spt ~ X2_febmay, data = master_data)

# Starry Flounder Bay Study + 10 ~ March - June X2 
stf_lm <- lm(log_stf ~ X2_marjun + step_1987, data = master_data)

# Bay Shrimp Bay Study ~ March - May X2
bsh_lm <- lm(log_bsh ~ X2_marmay + step_1987, data = master_data)
summary(bsh_lm)

lfs_lm$model
