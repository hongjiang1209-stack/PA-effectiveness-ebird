# 9.Complementarity.R
# 
# Purpose:
# Calculate protected-area complementarity based on the species observed in matched
# PA checklists and species-level observational scarcity weights.
#
# Input:
# Species weight table, combined matching results for 2010 and 2020, WDPAID lists,
# and matched eBird observation records.
#
# Output:
# RDS and CSV files containing WDPAID, number of species, and complementarity values
# for each protected area.
#
# Software:
# R 4.3.2
#
# Notes:
# Complementarity was calculated as the sum of species-level observational scarcity
# weights for species observed in matched PA checklists.

rm(list=ls())

species_weight<-readRDS("J:/us_ebird/Data_processing/7-1.Checklist_richw/species_weight.rds")
species_weight$Freq<-as.numeric(species_weight$Freq)


WDPA_2010<-readRDS("J:/us_ebird/Data_processing/Matching_2010/Matching_results_combined.rds")
WDPA_2010<-WDPA_2010[,1:2]
WDPA_2020<-readRDS("J:/us_ebird/Data_processing/Matching_2020/Matching_results_combined.rds")
WDPA_2020<-WDPA_2020[,1:2]
WDPA_chlist<-rbind(WDPA_2010,WDPA_2020)
WDPA_chlist<-as.data.frame((WDPA_chlist))
rm(WDPA_2010)
rm(WDPA_2020)


##### WDPA list
List_2010<-readRDS("J:/us_ebird/Data_processing/Matching_2010/WDPAID_nofchlist.rds")
List_2020<-readRDS("J:/us_ebird/Data_processing/Matching_2020/WDPAID_nofchlist.rds")
List_2010<-List_2010[,1:2]
List_2020<-List_2020[,1:2]
List_2010<-subset(List_2010,List_2010$nofchlist_PA>0)
List_2020<-subset(List_2020,List_2020$nofchlist_PA>0)
List<-rbind(List_2010,List_2020)
List<-as.data.frame(unique(List$WDPAID))
colnames(List)<-"WDPAID"

result<-data.frame(matrix(NA,nrow(List),3))
colnames(result)<-c("WDPAID","Nofspecies","freq_min")



records<-readRDS("J:/us_ebird/Data_processing/8.Complementarity/Matched_records_all.rds")
records<-records[,-(1:4)]
records<-records[,-(2:27)]

for(i in 1:nrow(List)){
  ID<-List$WDPAID[i]
  result$WDPAID[i]<-ID
  checklist_thisPA<-subset(WDPA_chlist,WDPA_chlist$WDPAID==ID)
  checklist_thisPA<-checklist_thisPA$checklist_PA
  records_thisPA<-subset(records,records$checklist %in% checklist_thisPA)
  species_thisPA<-unique(records_thisPA$scientific_name)
  result$Nofspecies[i]<-length(species_thisPA)
  Freq_thisPA<-subset(species_weight,species_weight$scientific_name %in% species_thisPA)
  result$freq_min[i]<-min(Freq_thisPA$Freq)
  print(paste0(i," in ",nrow(List)," FINISHED"))
}

saveRDS(result,"J:/us_ebird/Data_processing/8.Complementarity/PAs_minFreq.rds")



Complementarity<-data.frame(matrix(NA,nrow(result),3))
colnames(Complementarity)<-c("WDPAID","Nofspecies","complementarity")

for(i in 1:nrow(result)){
  ID<-List$WDPAID[i]
  Complementarity$WDPAID[i]<-ID
  checklist_thisPA<-subset(WDPA_chlist,WDPA_chlist$WDPAID==ID)
  checklist_thisPA<-checklist_thisPA$checklist_PA
  records_thisPA<-subset(records,records$checklist %in% checklist_thisPA)
  species_thisPA<-unique(records_thisPA$scientific_name)
  Freq_thisPA<-subset(species_weight,species_weight$scientific_name %in% species_thisPA)
  Complementarity$Nofspecies[i]<-length(species_thisPA)
  Complementarity$complementarity[i]<-sum(Freq_thisPA$weighte5)
  print(paste0(i," in ",nrow(result)," FINISHED"))
}

saveRDS(Complementarity,"J:/us_ebird/Data_processing/8.Complementarity/Complementarity.rds")
write.csv(Complementarity,"J:/us_ebird/Data_processing/8.Complementarity/Complementarity.csv")
