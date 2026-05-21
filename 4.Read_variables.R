# 4.Read_variables.R
# Purpose:
# Extract environmental, topographic, and socio-economic covariates for filtered eBird checklists,
# remove checklists with missing covariate values, and calculate z-score standardized covariates
# within each analysis period.
#
# Input:
# Filtered checklist-level RDS files for 2009-2011 and 2019-2021, and raster layers for NPP,
# distance to waterbody, elevation, TPI, construction land, and human population.
#
# Output:
# Checklist-level RDS files with raw and standardized covariates for 2010 and 2020.
#
#
# Notes: Script run on the University of York's HPC, Viking2.


# config ------------------------------------------------------------------

# Set working directory
directory <- '.' 


################################
###### array with realms #######
### 1 - $SLURM_ARRAY_TASK_ID ###
################################
cargs <- as.numeric(commandArgs(trailingOnly = TRUE))
array_index <- cargs[1]


time_list<-c("2010","2020")
time<-time_list[array_index]

# Load packages
print("Loading libraries")
library(raster);library(terra)
print('libraries loaded')




# Set location for outputted data to be written to 
folder_save <- file.path(directory,paste0("variables_",time))

if ( ! dir.exists(folder_save)){
  dir.create(folder_save, showWarnings = FALSE)
}

if(array_index==1){
  checklist<-readRDS("US2009_2011_chlist.rds")
}else{
  checklist<-readRDS("US2019_2021_chlist.rds")
}


print("number of checklist")
print(nrow(checklist))


#variables
NPP<-rast(paste0("MATCH_NPP_",time,".tif"))
DisToWater<-rast(paste0("MATCH_WATER_",time,".tif"))
Elevation<-rast("DEM_US.tif");
TPI<-rast("TPI_US.tif");
URBAN<-rast(paste0("MATCH_Urban_",time,".tif"))
POP<-rast(paste0("MATCH_POP_",time,".tif"))


#add columns
checklist[ , 'NPP']=NA
checklist[ , 'Elevation']=NA
checklist[ , 'URBAN']=NA
checklist[ , 'POP']=NA
checklist[ , 'TPI']=NA
checklist[ , 'Distowater']=NA


for(i in 1:nrow(checklist)){
  point<-cbind(checklist$lon[i],checklist$lat[i]);#lon,lat
  checklist$NPP[i]<-extract(NPP,point);
  checklist$Elevation[i]<-extract(Elevation,point);
  checklist$TPI[i]<-extract(TPI,point);
  checklist$Distowater[i]<-extract(DisToWater,point);
  checklist$URBAN[i]<-extract(URBAN,point);
  checklist$POP[i]<-extract(POP,point);
  print(paste0(i,"/",nrow(checklist)))
}

print("--------Variables read--------")

ch<-as.data.frame(checklist)
# ch <- ch[complete.cases(ch[, c("NPP", "Elevation", "TPI", "Distowater", "URBAN", "POP")]), ]
ch <- ch[complete.cases(ch[, c("NPP", "Elevation", "TPI", "Distowater")]), ]
saveRDS(ch,paste0(folder_save,"/Checklist_variables_noNA.rds_",time,".rds"))
checklist<-ch
checklist$NPP<-as.numeric(checklist$NPP)
checklist$Elevation<-as.numeric(checklist$Elevation)
checklist$TPI<-as.numeric(checklist$TPI)
checklist$Distowater<-as.numeric(checklist$Distowater)
checklist$URBAN<-as.numeric(checklist$URBAN)
checklist$POP<-as.numeric(checklist$POP)

#add z-score columns
checklist[ , 'z_lat']=NA
checklist[ , 'z_lon']=NA
checklist[ , 'z_NPP']=NA
checklist[ , 'z_Elevation']=NA
checklist[ , 'z_TPI']=NA
checklist[ , 'z_Distowater']=NA
checklist[ , 'z_URBAN']=NA
checklist[ , 'z_POP']=NA


#calculate mean and sd, for z-score normalization
mean_lat<-mean(checklist$lat);                 sd_lat<-sd(checklist$lat)
mean_lon<-mean(checklist$lon);                 sd_lon<-sd(checklist$lon)
mean_NPP<-mean(checklist$NPP);                 sd_NPP<-sd(checklist$NPP)
mean_Elevation<-mean(checklist$Elevation);     sd_Elevation<-sd(checklist$Elevation)
mean_TPI<-mean(checklist$TPI);                 sd_TPI<-sd(checklist$TPI)
mean_Distowater<-mean(checklist$Distowater);   sd_Distowater<-sd(checklist$Distowater)
# mean_URBAN<-mean(checklist$URBAN);             sd_URBAN<-sd(checklist$URBAN)
# mean_POP<-mean(checklist$POP);                 sd_POP<-sd(checklist$POP)

#z-score normalization
checklist$z_lat<-(checklist$lat-mean_lat)/sd_lat
checklist$z_lon<-(checklist$lon-mean_lon)/sd_lon
checklist$z_NPP<-(checklist$NPP-mean_NPP)/sd_NPP
checklist$z_Elevation<-(checklist$Elevation-mean_Elevation)/sd_Elevation
checklist$z_TPI<-(checklist$TPI-mean_TPI)/sd_TPI
checklist$z_Distowater<-(checklist$Distowater-mean_Distowater)/sd_Distowater
# checklist$z_URBAN<-(checklist$URBAN-mean_URBAN)/sd_URBAN
# checklist$z_POP<-(checklist$POP-mean_POP)/sd_POP

print("--------Z-score normalization finished--------")

#save
rownames(checklist)<-1:nrow(checklist)
saveRDS(checklist,paste0(folder_save,"/Checklist_variables_",time,".rds"))


