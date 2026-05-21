# 1.readEBD.R
# Purpose:
# Read raw eBird Basic Dataset (EBD) text files, retain approved complete checklists,
# remove unused fields, and save processed files as RDS objects.
#
# Input:
# Raw EBD text files downloaded from eBird.
#
# Output:
# RDS files containing filtered eBird observation records.
#
# Software:
# R 4.3.2
#
# Notes:
# Raw EBD files are not included in this repository due to data size and eBird data-use policies.
#
# This script was adapted from previously published eBird preprocessing workflows
# described in [Cazalis et al., 2020], with modifications for this study.
# Cazalis V, Princé K, Mihoub J, et al. 2020. Effectiveness of protected areas in conserving tropical forest birds. Nature Communications, 11, 4461.




rm(list = ls())

library(auk)

files.regions<-list.files("J:/ebird/txt", recursive=TRUE)

for (i in 1:length(files.regions)){
  
  eb<-read_ebd(paste0("J:/ebird/txt/",files.regions[i]), rollup=T, unique=F)
  print("read")
  # Checklist id (not the same than in script 2 because here we did not use the auk_unique function as we do not want to merge checklists of different observers)
  eb$checklist<-eb$sampling_event_identifier
  eb$checklist_id<-eb$sampling_event_identifier<-NULL
  print("eb$checklist_id")
  ### Reduce the dataset
  # Remove useless columns
  eb[,c("state_province", "subnational1_code","subspecies_scientific_name","state","state_code","global_unique_identifier", "last_edited_date", "subspecies_common_name", "breeding_bird_atlas_category", "effort_area_ha", "county", "subnational2_code", "iba_code", "bcr_code", "usfws_code", "atlas_block", "locality", "locality_id", "locality_type", "first_name", "last_name", "has_media", "reviewed", "x")]<-NULL
  # Remove disapproved observations
  eb<-subset(eb, eb$approved==T)
  print("approved")
  # Remove checklists that did not report all species detected
  eb<-subset(eb, eb$all_species_reported==T)
  print("all_species_reported")
  # Remove old observations
  eb$year<-as.numeric(format(eb$observation_date, "%Y"))
  print("year")
  saveRDS(eb, file=paste0("J:/ebird/rds/US_ebird_Export200123_", substr(files.regions[i],5,6), ".rds", sep=""))
  
  cat(i)
  print(i)
  
}
