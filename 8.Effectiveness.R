# 8.Effectiveness.R
# Purpose:
# Calculate local and spillover conservation effectiveness metrics for protected areas
# based on species richness estimates in PA, spillover, and control-area checklist groups.
#
# Input:
# Richness-duration estimation results for 2010 and 2020 across different richness metrics,
# including overall richness (richness), observational scarcity-weighted richness (richness weighted),
# (near-)threatened species richness (richness IUCN), and spatial rarity-weighted richness (richness rarity).
#
# Output:
# RDS files containing spatial, temporal, and spatio-temporal effectiveness estimates
# at local and spillover scales for each protected area.
#
# Software:
# R 4.3.2
#
# Notes:
# Effectiveness was calculated using both observed average richness and richness standardized
# to a 3-hour observation duration. The 3-hour standardized estimates were used in the main analyses.


rm(list=ls())
time_list<-c("2010","2020")

#   PA with matching result
WDPAID_2010<-readRDS("J:/us_ebird/Data_processing/Matching_2010/WDPAID_nofchlist.rds")
WDPAID_2020<-readRDS("J:/us_ebird/Data_processing/Matching_2020/WDPAID_nofchlist.rds")

WDPAID_2010<-subset(WDPAID_2010,WDPAID_2010$nofchlist_PA!=0)
WDPAID_2020<-subset(WDPAID_2020,WDPAID_2020$nofchlist_PA!=0)
WDPAID<-rbind(WDPAID_2010,WDPAID_2020)
WDPAID<-unique(WDPAID$WDPAID)  ## all the WDPAID


############################################################################
############################################################################
##                              Start richness                            ##
############################################################################
############################################################################

####################  start rich_avg  ######################
setwd("J:/us_ebird/Data_processing/7.Rich_duration")
WDPA_eff<-as.data.frame(WDPAID)
#richness
WDPA_eff[,"R_PA2010"]<-NA; WDPA_eff[,"R_PA2020"]<-NA; WDPA_eff[,"R_buffer2010"]<-NA
WDPA_eff[,"R_buffer2020"]<-NA; WDPA_eff[,"R_out2010"]<-NA; WDPA_eff[,"R_out2020"]<-NA

for(j in 1:2){ # read i
  time<-time_list[j] #time="2010" / "2020"
  rich<-readRDS(paste0("richness_duration",time,"_PBO_noW.rds"))
  
  rich_PA<-subset(rich,group=="PA")
  rich_buffer<-subset(rich,group=="buffer")
  rich_out<-subset(rich,group=="out")
  
  for(i in 1:nrow(WDPA_eff)){
    id<-WDPA_eff$WDPAID[i]
    
    ####################  start PA  ######################
    if(id %in% rich_PA$WDPAID){
      if(time=="2010"){
        WDPA_eff$R_PA2010[i]<-rich_PA$R_avg[which(rich_PA$WDPAID==id)]
      }else{
        WDPA_eff$R_PA2020[i]<-rich_PA$R_avg[which(rich_PA$WDPAID==id)]
      }
    }
    ####################  finish PA  ######################
    
    ####################  start buffer  ######################
    if(id %in% rich_buffer$WDPAID){
      if(time=="2010"){
        WDPA_eff$R_buffer2010[i]<-rich_buffer$R_avg[which(rich_buffer$WDPAID==id)]
      }else{
        WDPA_eff$R_buffer2020[i]<-rich_buffer$R_avg[which(rich_buffer$WDPAID==id)]
      }
    }
    ####################  finish buffer  ######################
    
    ####################  start out  ######################
    if(id %in% rich_out$WDPAID){
      if(time=="2010"){
        WDPA_eff$R_out2010[i]<-rich_out$R_avg[which(rich_out$WDPAID==id)]
      }else{
        WDPA_eff$R_out2020[i]<-rich_out$R_avg[which(rich_out$WDPAID==id)]
      }
    }
    ####################  finish out  ######################
    print(paste0(i," in ",nrow(WDPA_eff)," of rich_avg ",time," finished"))
  }
}


### calculate effectiveness of richness
WDPA_eff[,"Eff_R_S10_local"]<-(WDPA_eff$R_PA2010/WDPA_eff$R_out2010)-1
WDPA_eff[,"Eff_R_S20_local"]<-(WDPA_eff$R_PA2020/WDPA_eff$R_out2020)-1
WDPA_eff[,"Eff_R_S10_spill"]<-(WDPA_eff$R_buffer2010/WDPA_eff$R_out2010)-1
WDPA_eff[,"Eff_R_S20_spill"]<-(WDPA_eff$R_buffer2020/WDPA_eff$R_out2020)-1
WDPA_eff[,"Eff_R_T_local"]<-(WDPA_eff$R_PA2020/WDPA_eff$R_PA2010)-1
WDPA_eff[,"Eff_R_T_spill"]<-(WDPA_eff$R_buffer2020/WDPA_eff$R_buffer2010)-1
WDPA_eff[,"Eff_R_ST_local"]<-((WDPA_eff$R_PA2020/WDPA_eff$R_PA2010)/(WDPA_eff$R_out2020/WDPA_eff$R_out2010))-1
WDPA_eff[,"Eff_R_ST_spill"]<-((WDPA_eff$R_buffer2020/WDPA_eff$R_buffer2010)/(WDPA_eff$R_out2020/WDPA_eff$R_out2010))-1


saveRDS(WDPA_eff,"J:/us_ebird/Data_processing/9.Effectiveness/Eff_rich_avg.rds")
# WDPA_eff[is.na(WDPA_eff)]<-99999
# write.csv(WDPA_eff,"J:/us_ebird/Data_processing/9.Effectiveness/Eff_rich_avg.csv")
####################  finish rich_avg  ######################



####################  start rich_3H  ######################
WDPA_eff<-as.data.frame(WDPAID)
#richness
WDPA_eff[,"R_PA2010"]<-NA; WDPA_eff[,"R_PA2020"]<-NA; WDPA_eff[,"R_buffer2010"]<-NA
WDPA_eff[,"R_buffer2020"]<-NA; WDPA_eff[,"R_out2010"]<-NA; WDPA_eff[,"R_out2020"]<-NA

for(j in 1:2){ # read i
  time<-time_list[j] #time="2010" / "2020"
  rich<-readRDS(paste0("richness_duration",time,"_PBO_noW.rds"))
  
  rich_PA<-subset(rich,group=="PA")
  rich_buffer<-subset(rich,group=="buffer")
  rich_out<-subset(rich,group=="out")
  
  for(i in 1:nrow(WDPA_eff)){
    id<-WDPA_eff$WDPAID[i]
    
    ####################  start PA  ######################
    if(id %in% rich_PA$WDPAID){
      if(time=="2010"){
        WDPA_eff$R_PA2010[i]<-rich_PA$Estim_3H[which(rich_PA$WDPAID==id)]
      }else{
        WDPA_eff$R_PA2020[i]<-rich_PA$Estim_3H[which(rich_PA$WDPAID==id)]
      }
    }
    ####################  finish PA  ######################
    
    ####################  start buffer  ######################
    if(id %in% rich_buffer$WDPAID){
      if(time=="2010"){
        WDPA_eff$R_buffer2010[i]<-rich_buffer$Estim_3H[which(rich_buffer$WDPAID==id)]
      }else{
        WDPA_eff$R_buffer2020[i]<-rich_buffer$Estim_3H[which(rich_buffer$WDPAID==id)]
      }
    }
    ####################  finish buffer  ######################
    
    ####################  start out  ######################
    if(id %in% rich_out$WDPAID){
      if(time=="2010"){
        WDPA_eff$R_out2010[i]<-rich_out$Estim_3H[which(rich_out$WDPAID==id)]
      }else{
        WDPA_eff$R_out2020[i]<-rich_out$Estim_3H[which(rich_out$WDPAID==id)]
      }
    }
    ####################  finish out  ######################
    print(paste0(i," in ",nrow(WDPA_eff)," of rich_3H ",time," finished"))
  }
}

### calculate effectiveness of richness
WDPA_eff[,"Eff_R_S10_local"]<-(WDPA_eff$R_PA2010/WDPA_eff$R_out2010)-1
WDPA_eff[,"Eff_R_S20_local"]<-(WDPA_eff$R_PA2020/WDPA_eff$R_out2020)-1
WDPA_eff[,"Eff_R_S10_spill"]<-(WDPA_eff$R_buffer2010/WDPA_eff$R_out2010)-1
WDPA_eff[,"Eff_R_S20_spill"]<-(WDPA_eff$R_buffer2020/WDPA_eff$R_out2020)-1
WDPA_eff[,"Eff_R_T_local"]<-(WDPA_eff$R_PA2020/WDPA_eff$R_PA2010)-1
WDPA_eff[,"Eff_R_T_spill"]<-(WDPA_eff$R_buffer2020/WDPA_eff$R_buffer2010)-1
WDPA_eff[,"Eff_R_ST_local"]<-((WDPA_eff$R_PA2020/WDPA_eff$R_PA2010)/(WDPA_eff$R_out2020/WDPA_eff$R_out2010))-1
WDPA_eff[,"Eff_R_ST_spill"]<-((WDPA_eff$R_buffer2020/WDPA_eff$R_buffer2010)/(WDPA_eff$R_out2020/WDPA_eff$R_out2010))-1

saveRDS(WDPA_eff,"J:/us_ebird/Data_processing/9.Effectiveness/Eff_rich_3H.rds")
# WDPA_eff[is.na(WDPA_eff)]<-99999
# write.csv(WDPA_eff,"J:/us_ebird/Data_processing/9.Effectiveness/Eff_rich_3H.csv")
####################  finish rich_3H  ######################







############################################################################
############################################################################
##                        Start richness weighted                         ##
############################################################################
############################################################################

####################  start rich_avg_W ######################
setwd("J:/us_ebird/Data_processing/7.Rich_duration")
WDPA_eff<-as.data.frame(WDPAID)
#richness
WDPA_eff[,"R_PA2010"]<-NA; WDPA_eff[,"R_PA2020"]<-NA; WDPA_eff[,"R_buffer2010"]<-NA
WDPA_eff[,"R_buffer2020"]<-NA; WDPA_eff[,"R_out2010"]<-NA; WDPA_eff[,"R_out2020"]<-NA

for(j in 1:2){ # read i
  time<-time_list[j] #time="2010" / "2020"
  rich<-readRDS(paste0("richness_duration",time,"_PBO_W.rds"))
  
  rich_PA<-subset(rich,group=="PA")
  rich_buffer<-subset(rich,group=="buffer")
  rich_out<-subset(rich,group=="out")
  
  for(i in 1:nrow(WDPA_eff)){
    id<-WDPA_eff$WDPAID[i]
    
    ####################  start PA  ######################
    if(id %in% rich_PA$WDPAID){
      if(time=="2010"){
        WDPA_eff$R_PA2010[i]<-rich_PA$R_avg[which(rich_PA$WDPAID==id)]
      }else{
        WDPA_eff$R_PA2020[i]<-rich_PA$R_avg[which(rich_PA$WDPAID==id)]
      }
    }
    ####################  finish PA  ######################
    
    ####################  start buffer  ######################
    if(id %in% rich_buffer$WDPAID){
      if(time=="2010"){
        WDPA_eff$R_buffer2010[i]<-rich_buffer$R_avg[which(rich_buffer$WDPAID==id)]
      }else{
        WDPA_eff$R_buffer2020[i]<-rich_buffer$R_avg[which(rich_buffer$WDPAID==id)]
      }
    }
    ####################  finish buffer  ######################
    
    ####################  start out  ######################
    if(id %in% rich_out$WDPAID){
      if(time=="2010"){
        WDPA_eff$R_out2010[i]<-rich_out$R_avg[which(rich_out$WDPAID==id)]
      }else{
        WDPA_eff$R_out2020[i]<-rich_out$R_avg[which(rich_out$WDPAID==id)]
      }
    }
    ####################  finish out  ######################
    print(paste0(i," in ",nrow(WDPA_eff)," of rich_avg_W ",time," finished"))
  }
}

### calculate effectiveness of richness
WDPA_eff[,"Eff_R_S10_local"]<-(WDPA_eff$R_PA2010/WDPA_eff$R_out2010)-1
WDPA_eff[,"Eff_R_S20_local"]<-(WDPA_eff$R_PA2020/WDPA_eff$R_out2020)-1
WDPA_eff[,"Eff_R_S10_spill"]<-(WDPA_eff$R_buffer2010/WDPA_eff$R_out2010)-1
WDPA_eff[,"Eff_R_S20_spill"]<-(WDPA_eff$R_buffer2020/WDPA_eff$R_out2020)-1
WDPA_eff[,"Eff_R_T_local"]<-(WDPA_eff$R_PA2020/WDPA_eff$R_PA2010)-1
WDPA_eff[,"Eff_R_T_spill"]<-(WDPA_eff$R_buffer2020/WDPA_eff$R_buffer2010)-1
WDPA_eff[,"Eff_R_ST_local"]<-((WDPA_eff$R_PA2020/WDPA_eff$R_PA2010)/(WDPA_eff$R_out2020/WDPA_eff$R_out2010))-1
WDPA_eff[,"Eff_R_ST_spill"]<-((WDPA_eff$R_buffer2020/WDPA_eff$R_buffer2010)/(WDPA_eff$R_out2020/WDPA_eff$R_out2010))-1

saveRDS(WDPA_eff,"J:/us_ebird/Data_processing/9.Effectiveness/Eff_rich_avg_W.rds")

# WDPA_eff[is.na(WDPA_eff)]<-99999
# write.csv(WDPA_eff,"J:/us_ebird/Data_processing/9.Effectiveness/Eff_rich_avg_W.csv")
####################  finish rich_avg_W ######################



####################  start rich_3H_W  ######################
WDPA_eff<-as.data.frame(WDPAID)
#richness
WDPA_eff[,"R_PA2010"]<-NA; WDPA_eff[,"R_PA2020"]<-NA; WDPA_eff[,"R_buffer2010"]<-NA
WDPA_eff[,"R_buffer2020"]<-NA; WDPA_eff[,"R_out2010"]<-NA; WDPA_eff[,"R_out2020"]<-NA
for(j in 1:2){ # read i
  time<-time_list[j] #time="2010" / "2020"
  rich<-readRDS(paste0("richness_duration",time,"_PBO_W.rds"))
  
  rich_PA<-subset(rich,group=="PA")
  rich_buffer<-subset(rich,group=="buffer")
  rich_out<-subset(rich,group=="out")
  
  for(i in 1:nrow(WDPA_eff)){
    id<-WDPA_eff$WDPAID[i]
    
    ####################  start PA  ######################
    if(id %in% rich_PA$WDPAID){
      if(time=="2010"){
        WDPA_eff$R_PA2010[i]<-rich_PA$Estim_3H[which(rich_PA$WDPAID==id)]
      }else{
        WDPA_eff$R_PA2020[i]<-rich_PA$Estim_3H[which(rich_PA$WDPAID==id)]
      }
    }
    ####################  finish PA  ######################
    
    ####################  start buffer  ######################
    if(id %in% rich_buffer$WDPAID){
      if(time=="2010"){
        WDPA_eff$R_buffer2010[i]<-rich_buffer$Estim_3H[which(rich_buffer$WDPAID==id)]
      }else{
        WDPA_eff$R_buffer2020[i]<-rich_buffer$Estim_3H[which(rich_buffer$WDPAID==id)]
      }
    }
    ####################  finish buffer  ######################
    
    ####################  start out  ######################
    if(id %in% rich_out$WDPAID){
      if(time=="2010"){
        WDPA_eff$R_out2010[i]<-rich_out$Estim_3H[which(rich_out$WDPAID==id)]
      }else{
        WDPA_eff$R_out2020[i]<-rich_out$Estim_3H[which(rich_out$WDPAID==id)]
      }
    }
    ####################  finish out  ######################
    print(paste0(i," in ",nrow(WDPA_eff)," of rich_3H_W ",time," finished"))
  }
}

### calculate effectiveness of richness
WDPA_eff[,"Eff_R_S10_local"]<-(WDPA_eff$R_PA2010/WDPA_eff$R_out2010)-1
WDPA_eff[,"Eff_R_S20_local"]<-(WDPA_eff$R_PA2020/WDPA_eff$R_out2020)-1
WDPA_eff[,"Eff_R_S10_spill"]<-(WDPA_eff$R_buffer2010/WDPA_eff$R_out2010)-1
WDPA_eff[,"Eff_R_S20_spill"]<-(WDPA_eff$R_buffer2020/WDPA_eff$R_out2020)-1
WDPA_eff[,"Eff_R_T_local"]<-(WDPA_eff$R_PA2020/WDPA_eff$R_PA2010)-1
WDPA_eff[,"Eff_R_T_spill"]<-(WDPA_eff$R_buffer2020/WDPA_eff$R_buffer2010)-1
WDPA_eff[,"Eff_R_ST_local"]<-((WDPA_eff$R_PA2020/WDPA_eff$R_PA2010)/(WDPA_eff$R_out2020/WDPA_eff$R_out2010))-1
WDPA_eff[,"Eff_R_ST_spill"]<-((WDPA_eff$R_buffer2020/WDPA_eff$R_buffer2010)/(WDPA_eff$R_out2020/WDPA_eff$R_out2010))-1

saveRDS(WDPA_eff,"J:/us_ebird/Data_processing/9.Effectiveness/Eff_rich_3H_W.rds")
# WDPA_eff[is.na(WDPA_eff)]<-99999
# write.csv(WDPA_eff,"J:/us_ebird/Data_processing/9.Effectiveness/Eff_rich_3H_W.csv")
####################  finish rich_3H_W  ######################




############################################################################
############################################################################
##                          Start richness IUCN                           ##
############################################################################
############################################################################

####################  start rich_avg_IUCN ######################
setwd("J:/us_ebird/Data_processing/7.Rich_duration")
WDPA_eff<-as.data.frame(WDPAID)
#richness
WDPA_eff[,"R_PA2010"]<-NA; WDPA_eff[,"R_PA2020"]<-NA; WDPA_eff[,"R_buffer2010"]<-NA
WDPA_eff[,"R_buffer2020"]<-NA; WDPA_eff[,"R_out2010"]<-NA; WDPA_eff[,"R_out2020"]<-NA


for(j in 1:2){ # read i
  time<-time_list[j] #time="2010" / "2020"
  rich<-readRDS(paste0("richness_duration",time,"_PBO_IUCN.rds"))
  
  rich_PA<-subset(rich,group=="PA")
  rich_buffer<-subset(rich,group=="buffer")
  rich_out<-subset(rich,group=="out")
  
  for(i in 1:nrow(WDPA_eff)){
    id<-WDPA_eff$WDPAID[i]
    
    ####################  start PA  ######################
    if(id %in% rich_PA$WDPAID){
      if(time=="2010"){
        WDPA_eff$R_PA2010[i]<-rich_PA$R_avg[which(rich_PA$WDPAID==id)]
      }else{
        WDPA_eff$R_PA2020[i]<-rich_PA$R_avg[which(rich_PA$WDPAID==id)]
      }
    }
    ####################  finish PA  ######################
    
    ####################  start buffer  ######################
    if(id %in% rich_buffer$WDPAID){
      if(time=="2010"){
        WDPA_eff$R_buffer2010[i]<-rich_buffer$R_avg[which(rich_buffer$WDPAID==id)]
      }else{
        WDPA_eff$R_buffer2020[i]<-rich_buffer$R_avg[which(rich_buffer$WDPAID==id)]
      }
    }
    ####################  finish buffer  ######################
    
    ####################  start out  ######################
    if(id %in% rich_out$WDPAID){
      if(time=="2010"){
        WDPA_eff$R_out2010[i]<-rich_out$R_avg[which(rich_out$WDPAID==id)]
      }else{
        WDPA_eff$R_out2020[i]<-rich_out$R_avg[which(rich_out$WDPAID==id)]
      }
    }
    ####################  finish out  ######################
    print(paste0(i," in ",nrow(WDPA_eff)," of rich_avg_IUCN ",time," finished"))
  }
}

### calculate effectiveness of richness
WDPA_eff[,"Eff_R_S10_local"]<-(WDPA_eff$R_PA2010/WDPA_eff$R_out2010)-1
WDPA_eff[,"Eff_R_S20_local"]<-(WDPA_eff$R_PA2020/WDPA_eff$R_out2020)-1
WDPA_eff[,"Eff_R_S10_spill"]<-(WDPA_eff$R_buffer2010/WDPA_eff$R_out2010)-1
WDPA_eff[,"Eff_R_S20_spill"]<-(WDPA_eff$R_buffer2020/WDPA_eff$R_out2020)-1
WDPA_eff[,"Eff_R_T_local"]<-(WDPA_eff$R_PA2020/WDPA_eff$R_PA2010)-1
WDPA_eff[,"Eff_R_T_spill"]<-(WDPA_eff$R_buffer2020/WDPA_eff$R_buffer2010)-1
WDPA_eff[,"Eff_R_ST_local"]<-((WDPA_eff$R_PA2020/WDPA_eff$R_PA2010)/(WDPA_eff$R_out2020/WDPA_eff$R_out2010))-1
WDPA_eff[,"Eff_R_ST_spill"]<-((WDPA_eff$R_buffer2020/WDPA_eff$R_buffer2010)/(WDPA_eff$R_out2020/WDPA_eff$R_out2010))-1

saveRDS(WDPA_eff,"J:/us_ebird/Data_processing/9.Effectiveness/Eff_rich_avg_IUCN.rds")
# WDPA_eff[is.na(WDPA_eff)]<-99999
# WDPA_eff[sapply(WDPA_eff, is.infinite)] <- 0
# write.csv(WDPA_eff,"J:/us_ebird/Data_processing/9.Effectiveness/Eff_rich_avg_IUCN.csv")
####################  finish rich_avg_IUCN ######################



####################  start rich_3H_IUCN  ######################
WDPA_eff<-as.data.frame(WDPAID)
#richness
WDPA_eff[,"R_PA2010"]<-NA; WDPA_eff[,"R_PA2020"]<-NA; WDPA_eff[,"R_buffer2010"]<-NA
WDPA_eff[,"R_buffer2020"]<-NA; WDPA_eff[,"R_out2010"]<-NA; WDPA_eff[,"R_out2020"]<-NA
for(j in 1:2){ # read i
  time<-time_list[j] #time="2010" / "2020"
  rich<-readRDS(paste0("richness_duration",time,"_PBO_IUCN.rds"))
  
  rich_PA<-subset(rich,group=="PA")
  rich_buffer<-subset(rich,group=="buffer")
  rich_out<-subset(rich,group=="out")
  
  for(i in 1:nrow(WDPA_eff)){
    id<-WDPA_eff$WDPAID[i]
    
    ####################  start PA  ######################
    if(id %in% rich_PA$WDPAID){
      if(time=="2010"){
        WDPA_eff$R_PA2010[i]<-rich_PA$Estim_3H[which(rich_PA$WDPAID==id)]
      }else{
        WDPA_eff$R_PA2020[i]<-rich_PA$Estim_3H[which(rich_PA$WDPAID==id)]
      }
    }
    ####################  finish PA  ######################
    
    ####################  start buffer  ######################
    if(id %in% rich_buffer$WDPAID){
      if(time=="2010"){
        WDPA_eff$R_buffer2010[i]<-rich_buffer$Estim_3H[which(rich_buffer$WDPAID==id)]
      }else{
        WDPA_eff$R_buffer2020[i]<-rich_buffer$Estim_3H[which(rich_buffer$WDPAID==id)]
      }
    }
    ####################  finish buffer  ######################
    
    ####################  start out  ######################
    if(id %in% rich_out$WDPAID){
      if(time=="2010"){
        WDPA_eff$R_out2010[i]<-rich_out$Estim_3H[which(rich_out$WDPAID==id)]
      }else{
        WDPA_eff$R_out2020[i]<-rich_out$Estim_3H[which(rich_out$WDPAID==id)]
      }
    }
    ####################  finish out  ######################
    print(paste0(i," in ",nrow(WDPA_eff)," of rich_3H_IUCN ",time," finished"))
  }
}

### calculate effectiveness of richness
WDPA_eff[,"Eff_R_S10_local"]<-(WDPA_eff$R_PA2010/WDPA_eff$R_out2010)-1
WDPA_eff[,"Eff_R_S20_local"]<-(WDPA_eff$R_PA2020/WDPA_eff$R_out2020)-1
WDPA_eff[,"Eff_R_S10_spill"]<-(WDPA_eff$R_buffer2010/WDPA_eff$R_out2010)-1
WDPA_eff[,"Eff_R_S20_spill"]<-(WDPA_eff$R_buffer2020/WDPA_eff$R_out2020)-1
WDPA_eff[,"Eff_R_T_local"]<-(WDPA_eff$R_PA2020/WDPA_eff$R_PA2010)-1
WDPA_eff[,"Eff_R_T_spill"]<-(WDPA_eff$R_buffer2020/WDPA_eff$R_buffer2010)-1
WDPA_eff[,"Eff_R_ST_local"]<-((WDPA_eff$R_PA2020/WDPA_eff$R_PA2010)/(WDPA_eff$R_out2020/WDPA_eff$R_out2010))-1
WDPA_eff[,"Eff_R_ST_spill"]<-((WDPA_eff$R_buffer2020/WDPA_eff$R_buffer2010)/(WDPA_eff$R_out2020/WDPA_eff$R_out2010))-1

saveRDS(WDPA_eff,"J:/us_ebird/Data_processing/9.Effectiveness/Eff_rich_3H_IUCN.rds")
# WDPA_eff[is.na(WDPA_eff)]<-99999
# WDPA_eff[sapply(WDPA_eff, is.infinite)] <- 0
# write.csv(WDPA_eff,"J:/us_ebird/Data_processing/9.Effectiveness/Eff_rich_3H_IUCN.csv")
####################  finish rich_3H_IUCN  ######################









############################################################################
############################################################################
##                         Start richness rarity                          ##
############################################################################
############################################################################


####################  start rich_avg_rarity  ######################
setwd("J:/us_ebird/Data_processing/7.Rich_duration")
WDPA_eff<-as.data.frame(WDPAID)
#richness
WDPA_eff[,"R_PA2010"]<-NA; WDPA_eff[,"R_PA2020"]<-NA; WDPA_eff[,"R_buffer2010"]<-NA
WDPA_eff[,"R_buffer2020"]<-NA; WDPA_eff[,"R_out2010"]<-NA; WDPA_eff[,"R_out2020"]<-NA


for(j in 1:2){ # read i
  time<-time_list[j] #time="2010" / "2020"
  rich<-readRDS(paste0("richness_duration",time,"_PBO_rarity.rds"))
  
  rich_PA<-subset(rich,group=="PA")
  rich_buffer<-subset(rich,group=="buffer")
  rich_out<-subset(rich,group=="out")
  
  for(i in 1:nrow(WDPA_eff)){
    id<-WDPA_eff$WDPAID[i]
    
    ####################  start PA  ######################
    if(id %in% rich_PA$WDPAID){
      if(time=="2010"){
        WDPA_eff$R_PA2010[i]<-rich_PA$R_avg[which(rich_PA$WDPAID==id)]
      }else{
        WDPA_eff$R_PA2020[i]<-rich_PA$R_avg[which(rich_PA$WDPAID==id)]
      }
    }
    ####################  finish PA  ######################
    
    ####################  start buffer  ######################
    if(id %in% rich_buffer$WDPAID){
      if(time=="2010"){
        WDPA_eff$R_buffer2010[i]<-rich_buffer$R_avg[which(rich_buffer$WDPAID==id)]
      }else{
        WDPA_eff$R_buffer2020[i]<-rich_buffer$R_avg[which(rich_buffer$WDPAID==id)]
      }
    }
    ####################  finish buffer  ######################
    
    ####################  start out  ######################
    if(id %in% rich_out$WDPAID){
      if(time=="2010"){
        WDPA_eff$R_out2010[i]<-rich_out$R_avg[which(rich_out$WDPAID==id)]
      }else{
        WDPA_eff$R_out2020[i]<-rich_out$R_avg[which(rich_out$WDPAID==id)]
      }
    }
    ####################  finish out  ######################
    print(paste0(i," in ",nrow(WDPA_eff)," of rich_avg_rarity ",time," finished"))
  }
}


### calculate effectiveness of richness
WDPA_eff[,"Eff_R_S10_local"]<-(WDPA_eff$R_PA2010/WDPA_eff$R_out2010)-1
WDPA_eff[,"Eff_R_S20_local"]<-(WDPA_eff$R_PA2020/WDPA_eff$R_out2020)-1
WDPA_eff[,"Eff_R_S10_spill"]<-(WDPA_eff$R_buffer2010/WDPA_eff$R_out2010)-1
WDPA_eff[,"Eff_R_S20_spill"]<-(WDPA_eff$R_buffer2020/WDPA_eff$R_out2020)-1
WDPA_eff[,"Eff_R_T_local"]<-(WDPA_eff$R_PA2020/WDPA_eff$R_PA2010)-1
WDPA_eff[,"Eff_R_T_spill"]<-(WDPA_eff$R_buffer2020/WDPA_eff$R_buffer2010)-1
WDPA_eff[,"Eff_R_ST_local"]<-((WDPA_eff$R_PA2020/WDPA_eff$R_PA2010)/(WDPA_eff$R_out2020/WDPA_eff$R_out2010))-1
WDPA_eff[,"Eff_R_ST_spill"]<-((WDPA_eff$R_buffer2020/WDPA_eff$R_buffer2010)/(WDPA_eff$R_out2020/WDPA_eff$R_out2010))-1



saveRDS(WDPA_eff,"J:/us_ebird/Data_processing/9.Effectiveness/Eff_rich_avg_rarity.rds")
# WDPA_eff[is.na(WDPA_eff)]<-99999
# WDPA_eff[sapply(WDPA_eff, is.infinite)] <- 0
# write.csv(WDPA_eff,"J:/us_ebird/Data_processing/9.Effectiveness/Eff_rich_avg_rarity.csv")
####################  finish rich_avg_rarity  ######################




####################  start rich_3H_rarity  ######################
WDPA_eff<-as.data.frame(WDPAID)
#richness
WDPA_eff[,"R_PA2010"]<-NA; WDPA_eff[,"R_PA2020"]<-NA; WDPA_eff[,"R_buffer2010"]<-NA
WDPA_eff[,"R_buffer2020"]<-NA; WDPA_eff[,"R_out2010"]<-NA; WDPA_eff[,"R_out2020"]<-NA

for(j in 1:2){ # read i
  time<-time_list[j] #time="2010" / "2020"
  rich<-readRDS(paste0("richness_duration",time,"_PBO_rarity.rds"))
  
  rich_PA<-subset(rich,group=="PA")
  rich_buffer<-subset(rich,group=="buffer")
  rich_out<-subset(rich,group=="out")
  
  for(i in 1:nrow(WDPA_eff)){
    id<-WDPA_eff$WDPAID[i]
    
    ####################  start PA  ######################
    if(id %in% rich_PA$WDPAID){
      if(time=="2010"){
        WDPA_eff$R_PA2010[i]<-rich_PA$Estim_3H[which(rich_PA$WDPAID==id)]
      }else{
        WDPA_eff$R_PA2020[i]<-rich_PA$Estim_3H[which(rich_PA$WDPAID==id)]
      }
    }
    ####################  finish PA  ######################
    
    ####################  start buffer  ######################
    if(id %in% rich_buffer$WDPAID){
      if(time=="2010"){
        WDPA_eff$R_buffer2010[i]<-rich_buffer$Estim_3H[which(rich_buffer$WDPAID==id)]
      }else{
        WDPA_eff$R_buffer2020[i]<-rich_buffer$Estim_3H[which(rich_buffer$WDPAID==id)]
      }
    }
    ####################  finish buffer  ######################
    
    ####################  start out  ######################
    if(id %in% rich_out$WDPAID){
      if(time=="2010"){
        WDPA_eff$R_out2010[i]<-rich_out$Estim_3H[which(rich_out$WDPAID==id)]
      }else{
        WDPA_eff$R_out2020[i]<-rich_out$Estim_3H[which(rich_out$WDPAID==id)]
      }
    }
    ####################  finish out  ######################
    print(paste0(i," in ",nrow(WDPA_eff)," of rich_3H_rarity ",time," finished"))
  }
}

### calculate effectiveness of richness
WDPA_eff[,"Eff_R_S10_local"]<-(WDPA_eff$R_PA2010/WDPA_eff$R_out2010)-1
WDPA_eff[,"Eff_R_S20_local"]<-(WDPA_eff$R_PA2020/WDPA_eff$R_out2020)-1
WDPA_eff[,"Eff_R_S10_spill"]<-(WDPA_eff$R_buffer2010/WDPA_eff$R_out2010)-1
WDPA_eff[,"Eff_R_S20_spill"]<-(WDPA_eff$R_buffer2020/WDPA_eff$R_out2020)-1
WDPA_eff[,"Eff_R_T_local"]<-(WDPA_eff$R_PA2020/WDPA_eff$R_PA2010)-1
WDPA_eff[,"Eff_R_T_spill"]<-(WDPA_eff$R_buffer2020/WDPA_eff$R_buffer2010)-1
WDPA_eff[,"Eff_R_ST_local"]<-((WDPA_eff$R_PA2020/WDPA_eff$R_PA2010)/(WDPA_eff$R_out2020/WDPA_eff$R_out2010))-1
WDPA_eff[,"Eff_R_ST_spill"]<-((WDPA_eff$R_buffer2020/WDPA_eff$R_buffer2010)/(WDPA_eff$R_out2020/WDPA_eff$R_out2010))-1


saveRDS(WDPA_eff,"J:/us_ebird/Data_processing/9.Effectiveness/Eff_rich_3H_rarity.rds")
# WDPA_eff[is.na(WDPA_eff)]<-99999
# WDPA_eff[sapply(WDPA_eff, is.infinite)] <- 0
# write.csv(WDPA_eff,"J:/us_ebird/Data_processing/9.Effectiveness/Eff_rich_3H_rarity.csv")
####################  finish rich_3H_rarity  ######################


