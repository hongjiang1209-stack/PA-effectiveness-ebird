# 3.number_of_checklists_in_each_PA.R
# Purpose:
# Count the number of filtered eBird checklists located within each protected area
# for the 2010 and 2020 analysis periods.
#
# Input:
# Protected area polygon shpfile and filtered checklist-level RDS files.
#
# Output:
# CSV files containing WDPAID, protected area area, and the number of checklists
# within each protected area.
#
# Software:
# R 4.3.2
#
# Notes:
# File paths should be modified according to the local directory structure before reuse.


rm(list=ls())

library(sp); library(raster); library(sf); library(dplyr); library(stats)


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



Times<-c("2010","2020")
for(s in 1:2){
  time<-Times[s]
  
  setwd("J:/us_ebird/Data_processing/2.Nofchlist_inPA")
  Polygon_US<-st_read("J:/us_ebird/Data_processing/5.Matching/PA_us.shp")

  if(s==1){
  checklist<-readRDS("J:/us_ebird/Data_processing/1.Merge&filter/filter_100_10_15/US_2009_2011_chlist.rds")
  }else{
    checklist<-readRDS("J:/us_ebird/Data_processing/1.Merge&filter/filter_100_10_15/US_2019_2021_chlist.rds")
  }
  
  
  #an empty matrix saving WDPAID and number of checklists
  columns<- c("OBJECTID","WDPAID","NofChlist","area")
  result<- data.frame(matrix(nrow = nrow(Polygon_US), ncol = length(columns)))
  colnames(result)<-columns
  
  
  for(i in 1:nrow(Polygon_US)){
  #for(i in 1:5){
    poly<-Polygon_US[i,] #read the i th polygon
    num<-0;#number of checklist, from 0 for each PA
    shp<-Polygon_US[i,]$geometry[[1]];
    extent_poly<-extent_multipoly(shp);#return xmin, xmax, ymin, ymax
    #if no checklist in this extent
    a<-(checklist$lon>=extent_poly$xmin)&(checklist$lon<=extent_poly$xmax)
    sub_checklist<-subset(checklist,a)
    b<-(sub_checklist$lat>=extent_poly$ymin)&(sub_checklist$lat<=extent_poly$ymax)
    sub_checklist<-subset(sub_checklist,b)
    
    if(nrow(sub_checklist)!=0){
      for(j in 1:nrow(sub_checklist)){
        #subset checklist
        if(point.in.multipolygon(sub_checklist$lat[j],sub_checklist$lon[j],Polygon_US[i,])!=0){
          num<-num+1;
        }
      }
      }
    result[i,1]<-i
    result[i,2]<-Polygon_US[i,]$WDPAID
    result[i,3]<-num;
    result[i,4]<-Polygon_US[i,]$area
    print(paste(i,"/",nrow(Polygon_US)))
  }
  
  #output
  write.csv(result, paste0("number of checklist in each PA_",time,".csv"), row.names=FALSE)
  print(paste0(time," finished"))
}
