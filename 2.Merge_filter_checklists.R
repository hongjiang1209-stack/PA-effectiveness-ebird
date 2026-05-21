# 2.Merge_filter_checklists.R
# Purpose:
# Merge processed eBird records, generate checklist-level summaries, identify experienced observers,
# apply checklist filters, remove potential duplicate checklists, and export filtered records for
# the 2009-2011 and 2019-2021 analysis periods.
#
# Input:
# RDS files produced from raw eBird Basic Dataset records.
#
# Output:
# Filtered checklist-level and record-level RDS files for 2009-2011 and 2019-2021.
#
# Software:
# R 4.3.2
#
# Notes:
# File paths should be modified according to the local directory structure before reuse.


rm(list = ls())

time="2001_2023"

`%not in%` <- function (x, table) is.na(match(x, table, nomatch=NA_integer_)) 
#load the libraries
library(auk) ; library(plyr) ; library(readr) ; library(nlme);library(dplyr)
library(sf); library(caret);library(base);library(progress)

###################
### Merge Files ###
###################

#set the work path

# #####################  US-49 states #####################
setwd("J:/ebird/rds")
output_path="J:/ebird/checklist"
# #####################  US-49 states #####################


df<-readRDS(list.files()[1])
print(list.files()[1])
if(length(list.files())>1){
for(i in 2:length(list.files())){
  df_add<-readRDS(list.files()[i])
  df<-rbind(df,df_add)
  print(list.files()[i])
}
}


chlist<-ddply(df, .(df$checklist), function(x){data.frame(
  date=x$observation_date[1],
  time=x$time_observations_started[1],
  distance=x$effort_distance_km[1],
  n.observers=x$number_observers[1],
  duration=x$duration_minutes[1],
  observer=x$observer_id[1],
  rich=nrow(x),
  lon=x$longitude[1],
  lat=x$latitude[1],
  protocol=x$protocol_type[1],
  year=x$year[1],
  country=x$country[1]
)})

names(chlist)<-replace(names(chlist), names(chlist)=="df$checklist", "checklist")

chlist$day<-as.numeric(format(chlist$date, "%j"))

chlist$time.min<-sapply(strsplit(as.character(chlist$time),":"),
                        function(x) {
                          x <- as.numeric(x)
                          x[1]*60+x[2] })
saveRDS(chlist,paste0(output_path,"/chlist_mergedALL.rds"))

### Observer
observer.var<-ddply(df, .(df$observer_id), function(x){data.frame(Nb_obs=nrow(x), Nb_spc=nlevels(droplevels(as.factor(x$scientific_name))), Nb_checklist=nlevels(droplevels(as.factor(x$checklist))))})
names(observer.var)[1]<-"observer_id"
saveRDS(observer.var,paste0(output_path,"/observers_all.rds"))

######Processing filters

#Filter 1: Complete list(has been finished when processing the data of different states)

#Filter 2: experienced observers
# Remove observers with less than 100 species or less than 10 checklists and 15 species per checklist
obs.to.exclude<-subset(observer.var, observer.var$Nb_spc<100 | observer.var$Nb_checklist<10 | (observer.var$Nb_obs/observer.var$Nb_checklist)<15)
chlist2<-subset(chlist, chlist$observer %not in% obs.to.exclude$observer_id)
chlist2$observer<-droplevels(chlist2$observer)

print("chlist2")

#Filter 3: 30min<=duration<=600min
chlist3<-subset(chlist2,(chlist2$duration>=30)&(chlist2$duration<=600))
print("chlist3")

#Filter 4: distance<=10km
chlist4<-subset(chlist3,chlist3$distance<=10)
print("chlist4")

#Filter 5: Remove NAs
cat(paste(100*round(table(is.na(chlist4$duration))["TRUE"]/nrow(chlist4),2), " % of duration values are NA"))
chlist5<-chlist4[,c("checklist", "rich", "protocol", "duration", "n.observers", "time.min", "lon", "lat", "day", "observer", "year")]
chlist6<-chlist5[complete.cases(chlist5),]
chlist6$observer<-droplevels(as.factor(chlist6$observer))
print("chlist6")

#Filter 6: Moving potential duplicates (checklists made on the same day at the same place)
chlist6<-chlist6 %>%arrange(year,day,lon,lat,desc(rich))#arrange according to day, lat, lon, lat
chlist7<-distinct(chlist6,year,day,lon,lat,time.min,duration, .keep_all=TRUE)
print("chlist7")

#save checklists
saveRDS(chlist, file.path(output_path,paste0("US_",time,"_chlist.rds")))
saveRDS(chlist4, file.path(output_path,paste0("US_",time,"_chlist4.rds")))
saveRDS(chlist6, file.path(output_path,paste0("US_",time,"_chlist6.rds")))
saveRDS(chlist7, file.path(output_path,paste0("US_",time,"_chlist7.rds")))

#extract records based on chlist7
records<-subset(df,df$checklist %in% chlist7$checklist)

#output records
saveRDS(records, file.path(output_path,paste0("US_",time,"_records_filtered7.rds")))

##year 2009-2011; 2019-2021
chlist7$month <- as.numeric(format(chlist7$date, "%m"))

chlist7_breeding <- subset(chlist7, month %in% c(5, 6, 7))

chlist7_2010 <- subset(chlist7_breeding, year %in% c(2009, 2010, 2011))
chlist7_2020 <- subset(chlist7_breeding, year %in% c(2019, 2020, 2021))

records_2010 <- subset(df, checklist %in% chlist7_2010$checklist)
records_2020 <- subset(df, checklist %in% chlist7_2020$checklist)

saveRDS(chlist7_2010, file.path(output_path, "US_2009_2011_chlist.rds"))
saveRDS(chlist7_2020, file.path(output_path, "US_2019_2021_chlist.rds"))
saveRDS(records_2010, file.path(output_path, "US_2009_2011_records_filtered.rds"))
saveRDS(records_2020, file.path(output_path, "US_2019_2021_records_filtered.rds"))