# 5.Subset_checklists.R
# Purpose:
# Classify filtered eBird checklists into protected area (PA), spillover area (SA),
# and outside/control-area candidate groups based on protected area boundaries and
# 5-km buffer polygons.
#
# Input:
# Protected area polygon shpfile, 5-km buffer polygon shpfile, and checklist-level
# RDS files with extracted covariates.
#
# Output:
# RDS files containing checklists inside PAs, checklists within 5-km spillover buffers,
# and checklists outside both PAs and spillover buffers.
#
#
# Notes: Script run on the University of York's HPC, Viking2.



# config ------------------------------------------------------------------

# Set working directory
directory <- '.' 

# Load packages
print("Loading libraries")
###libraries
library(sp); library(raster); library(sf); library(dplyr)
print('libraries loaded')

################################
###### array with realms #######
### 1 - $SLURM_ARRAY_TASK_ID ###
################################
cargs <- as.numeric(commandArgs(trailingOnly = TRUE))
array_index <- cargs[1]

time_list<-c("2010","2020")
time<-time_list[array_index]

# Set location for outputted data to be written to 
folder_save <- file.path(directory,"results")

if ( ! dir.exists(folder_save)){
  dir.create(folder_save, showWarnings = FALSE)
}


#Define a function, to judge whether a point is inside a multi-element polygon
#input should be Polygon_US[i,]
point.in.multipolygon<-function(lat,lon,polygons){
  #library(sp);
  count<-0 #count of hole
  re<-0  #return value
  geo<-polygons$geometry[[1]]
  #one/multi_element
  for(i in 1:length(geo)){# number of elements
    #if the point is inside any single element, return 1
    if(!is.list(geo[[1]])){
      in_polygon<-point.in.polygon(lon,lat,geo[[i]][,1],geo[[i]][,2])
      if(in_polygon!=0){
        count<-count+1
      }
    }else{
      for(j in 1:length(geo[[i]])){#number of close polygon (including holes)
        #if inside number is double, not in this element; if inside number is single, in this element
        in_polygon<-point.in.polygon(lon,lat,geo[[i]][[j]][,1],geo[[i]][[j]][,2])
        if(in_polygon!=0){
          count<-count+1
        }
      }
    }
  }
  if(count%%2!=0){
    re<-1
  }
  return(re)
}



###############################################################
#define a function, get the extent of lat and lon of a multi-element polygon
#input should be Polygon_US[i,]$geometry[[1]]
extent_multipoly<-function(shp){
  if(length(shp)==1){#only one list
    if(!is.list(shp[[1]])){#a single element, no holes,[[1]] with the form of A*2
      Xmin<-extent(shp[[1]])@xmin
      Xmax<-extent(shp[[1]])@xmax
      Ymin<-extent(shp[[1]])@ymin
      Ymax<-extent(shp[[1]])@ymax 
    }else{#a single element
      #a single element, no holes，[[1]]with[[1]],A*2 inside
      Xmin<-extent(shp[[1]][[1]])@xmin
      Xmax<-extent(shp[[1]][[1]])@xmax
      Ymin<-extent(shp[[1]][[1]])@ymin
      Ymax<-extent(shp[[1]][[1]])@ymax
      if(length(shp[[1]])>1){##with hole
        for(i in 2:length(shp[[1]])){
          Xmin<-min(Xmin,extent(shp[[1]][[i]])@xmin);
          Xmax<-max(Xmax,extent(shp[[1]][[i]])@xmax);
          Ymin<-min(Ymin,extent(shp[[1]][[i]])@ymin);
          Ymax<-max(Ymax,extent(shp[[1]][[i]])@ymax);
        }
      }
    }
  }else{#multiple-element
    if(!is.list(shp[[1]])){
      Xmin<-extent(shp[[1]])@xmin
      Xmax<-extent(shp[[1]])@xmax
      Ymin<-extent(shp[[1]])@ymin
      Ymax<-extent(shp[[1]])@ymax 
    }else{
      Xmin<-extent(shp[[1]][[1]])@xmin
      Xmax<-extent(shp[[1]][[1]])@xmax
      Ymin<-extent(shp[[1]][[1]])@ymin
      Ymax<-extent(shp[[1]][[1]])@ymax
      if(length(shp[[1]])>1){
        for(i in 2:lengths(shp)[1]){
          Xmin<-min(Xmin,extent(shp[[1]][[i]])@xmin);
          Xmax<-max(Xmax,extent(shp[[1]][[i]])@xmax);
          Ymin<-min(Ymin,extent(shp[[1]][[i]])@ymin);
          Ymax<-max(Ymax,extent(shp[[1]][[i]])@ymax);
        }
      }
    }
    for(j in 2:length(shp)){
      if(!is.list(shp[[j]])){
        Xmin<-min(Xmin,extent(shp[[j]])@xmin);
        Xmax<-max(Xmax,extent(shp[[j]])@xmax);
        Ymin<-min(Ymin,extent(shp[[j]])@ymin);
        Ymax<-max(Ymax,extent(shp[[j]])@ymax);
      }else{
        for(m in 1:lengths(shp)[j]){
          Xmin<-min(Xmin,extent(shp[[j]][[m]])@xmin);
          Xmax<-max(Xmax,extent(shp[[j]][[m]])@xmax);
          Ymin<-min(Ymin,extent(shp[[j]][[m]])@ymin);
          Ymax<-max(Ymax,extent(shp[[j]][[m]])@ymax);
        }
      }
    }
  }
  
  out<-list(xmin=Xmin, xmax=Xmax, ymin=Ymin, ymax=Ymax)
  return(out)
}

print("Two functions defined")


#read data
Polygon_PA<-st_read("US_all.shp")
Polygon_buffer<-st_read("US_all_buffer5.shp")#merged and then seperated. each polygon has just one part

checklist<-readRDS(paste0("Checklist_variables_",time,".rds"))

checklist[ , 'PA'] = 0;
#for PA judge

for(i in 1:nrow(Polygon_PA)){
  poly<-Polygon_PA[i,]
  shp<-Polygon_PA[i,]$geometry[[1]]
  extent_PA<-extent_multipoly(shp);
  print("loop for checklist")
  for(j in 1:nrow(checklist)){
    if(checklist$PA[j]==0){
      if((checklist$lon[j]>=extent_PA$xmin)&(checklist$lon[j]<=extent_PA$xmax)){
        if((checklist$lat[j]>=extent_PA$ymin)&(checklist$lat[j]<=extent_PA$ymax)){
          checklist$PA[j]<-point.in.multipolygon(checklist$lat[j],checklist$lon[j],poly)
          #print(paste0(j," in ",i," th PA"))
        }
      }
    }
  }#PA==0 outside, else, inside
  
  #process report
  print(paste0(i,"/",nrow(Polygon_PA)," in PA"))
}


print("PA finished")
saveRDS(checklist,paste0("checklist_",time,"_judgePA.rds"))
chlist<-checklist
checklist_PA<-subset(chlist,chlist$PA==1)
checklist_notPA<-subset(chlist,chlist$PA==0)

colnames(checklist_notPA)[ncol(checklist_notPA)] <- "buffer"
for(i in 1:nrow(Polygon_buffer)){
  poly<-Polygon_buffer[i,]
  shp<-Polygon_buffer[i,]$geometry[[1]]
  extent_buffer<-extent_multipoly(shp)
  for(j in 1:nrow(checklist_notPA)){
    if(checklist_notPA$buffer[j]==0){
      if((checklist_notPA$lon[j]>=extent_buffer$xmin)&(checklist_notPA$lon[j]<=extent_buffer$xmax)){
        if((checklist_notPA$lat[j]>=extent_buffer$ymin)&(checklist_notPA$lat[j]<=extent_buffer$ymax)){
          checklist_notPA$buffer[j]<-point.in.multipolygon(checklist_notPA$lat[j],checklist_notPA$lon[j],poly)
          #print(paste0(j," in ",i," th buffer"))
        }
      }
    }
  }
  #process report
  print(paste0(i,"/",nrow(Polygon_buffer)," in buffer"))
}


checklist_buffer<-subset(checklist_notPA,checklist_notPA$buffer!=0)
checklist_out<-subset(checklist_notPA,checklist_notPA$buffer==0)
print("subset finished")

#save
saveRDS(checklist_PA, paste0(folder_save,"/checklist_PA_",time,".rds"))
saveRDS(checklist_buffer, paste0(folder_save,"/checklist_buffer_",time,".rds"))
saveRDS(checklist_out, paste0(folder_save,"/checklist_out_",time,".rds"))
print("F I N I S H E D")
