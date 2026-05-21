# 6.Matching_checklists.R
# Purpose:
# Match checklists within protected areas (PA), spillover areas (SA), and control-area
# candidates for each protected area using standardized geographic, environmental,
# and topographic covariates.
#
# Input:
# Checklist RDS files classified as PA, SA, and control-area candidates; protected area,
# 5-km buffer, and 5-50 km control-area shpfiles; and a WDPAID list.
#
# Output:
# CSV files containing matched checklist IDs for each protected area.
#
# Matching covariates:
# longitude, latitude, NPP, elevation, TPI, and distance to waterbody.
#
#
# Notes: Script run on the University of York's HPC, Viking2.



# config ------------------------------------------------------------------

# Set working directory
directory <- '.' 

# Load packages
print("Loading libraries")
###libraries
library(sp); library(raster); library(sf); library(dplyr); library(stats)
print('libraries loaded')

################################
###### array with realms #######
### 1 - $SLURM_ARRAY_TASK_ID ###
################################
cargs <- as.numeric(commandArgs(trailingOnly = TRUE))
array_index <- cargs[1]


# Set location for outputted data to be written to 
folder_save <- file.path(directory,paste0("matching_result"))

if ( ! dir.exists(folder_save)){
  dir.create(folder_save, showWarnings = FALSE)
}


#定义一个函数，判断是否在多要素多边形内
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
#定义一个函数，获取多要素多边形的经纬度范围（lat和lon的最大最小值范围）
#define a function, get the extent of lat and lon of a multi-element polygon
#input should be Polygon_US[i,]$geometry[[1]]
extent_multipoly<-function(shp){
  if(length(shp)==1){#只有一个list
    if(!is.list(shp[[1]])){#一个单独的元素，无孔,且一级[[1]]下就是A*2形式
      Xmin<-extent(shp[[1]])@xmin
      Xmax<-extent(shp[[1]])@xmax
      Ymin<-extent(shp[[1]])@ymin
      Ymax<-extent(shp[[1]])@ymax 
    }else{#一个单独的元素
      #如果一个元素，无孔，一级的[[1]]下为[[1]]，打开是A*2形式
      Xmin<-extent(shp[[1]][[1]])@xmin
      Xmax<-extent(shp[[1]][[1]])@xmax
      Ymin<-extent(shp[[1]][[1]])@ymin
      Ymax<-extent(shp[[1]][[1]])@ymax
      if(length(shp[[1]])>1){##如果该元素有孔
        for(i in 2:length(shp[[1]])){
          Xmin<-min(Xmin,extent(shp[[1]][[i]])@xmin);
          Xmax<-max(Xmax,extent(shp[[1]][[i]])@xmax);
          Ymin<-min(Ymin,extent(shp[[1]][[i]])@ymin);
          Ymax<-max(Ymax,extent(shp[[1]][[i]])@ymax);
        }
      }
    }
  }else{#多个要素
    if(!is.list(shp[[1]])){#第一个要素是无孔的,且[[1]]下面就是A*2
      Xmin<-extent(shp[[1]])@xmin
      Xmax<-extent(shp[[1]])@xmax
      Ymin<-extent(shp[[1]])@ymin
      Ymax<-extent(shp[[1]])@ymax 
    }else{#第一个要素
      #如果无孔，且[[1]]下面还有[[1]]，下面就是A*2
      Xmin<-extent(shp[[1]][[1]])@xmin
      Xmax<-extent(shp[[1]][[1]])@xmax
      Ymin<-extent(shp[[1]][[1]])@ymin
      Ymax<-extent(shp[[1]][[1]])@ymax
      if(length(shp[[1]])>1){##如果该元素有孔
        for(i in 2:lengths(shp)[1]){
          Xmin<-min(Xmin,extent(shp[[1]][[i]])@xmin);
          Xmax<-max(Xmax,extent(shp[[1]][[i]])@xmax);
          Ymin<-min(Ymin,extent(shp[[1]][[i]])@ymin);
          Ymax<-max(Ymax,extent(shp[[1]][[i]])@ymax);
        }
      }
    }
    #从第二个要素开始
    for(j in 2:length(shp)){
      if(!is.list(shp[[j]])){#第j个要素是无孔的,且[[1]]下面就是A*2
        Xmin<-min(Xmin,extent(shp[[j]])@xmin);
        Xmax<-max(Xmax,extent(shp[[j]])@xmax);
        Ymin<-min(Ymin,extent(shp[[j]])@ymin);
        Ymax<-max(Ymax,extent(shp[[j]])@ymax);
      }else{#第j个要素
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
##########################################################################################
##############读入需要输入PA和buffer的范围数据；有六个协变量字段的checklist###############
#####read the PA and buffer.shp; and checklists with six fields of matching variables#####
##########################################################################################

##Checklist_with six variables
checklist_PA<-readRDS("checklist_PA.rds")#读取read
checklist_buffer<-readRDS("checklist_buffer.rds")#读取read
checklist_out<-readRDS("checklist_out.rds")#读取read


##读取PA和buffer数据；WDPAID和buffer字段用于区分
##read the PAs and buffer data, "0"-PA,"1"-5km buffer,"2"-20km buffer
Polygon_PA<-st_read("PA_US.shp")
Polygon_buffer<-st_read("Buffer_US.shp")
Polygon_out<-st_read("Out_US.shp")
WDPAID<-read.csv("WDPAID.csv",header=FALSE)#all the WDPAID of the PAs, for the loop

####################  loop for each WDPAID  ###############################
## 1. read the extent of this PA and the corresponding buffer  ############
## 2. filter the checklists within the extent  ############################
## 3. get subsets of checklists from checklist_PA/buffer/out  #############
## 4. delete checklists outside the three polygons  #######################
## 5. matching, remove the checklists after matching, save ID  ############
## 5.1 calculate the mean values of variables #############################
## 5.2 checklist_PA>=/<checklist_buffer: if/else  #########################
## 5.3 arrange according to distance to mean values  ######################
## 5.4 calculate distance from the matching checklist  ####################
## 5.5 arrange checklists according to distance, desc(lat,lon)  ###########
## 5.6 (if) extract the first one, save to results, delete from subsets ###
## 5.6 (else) extract V_PA of the matching chlist, match out with PA  #####
## 5.7 (else) delete the chlist of PA and out  ############################
## 5.8 (esle) matching the remaining chlist_PA with out  ##################
###########################################################################


## 1. read PA/buffer/out of this WDPAID
ID<-WDPAID[array_index,1]
PA<-subset(Polygon_PA,Polygon_PA$WDPAID==ID)
buffer<-subset(Polygon_buffer,Polygon_buffer$WDPAID==ID)
out<-subset(Polygon_out,Polygon_out$WDPAID==ID)

## 2. filter the checklists within the extent
## 3. get subsets of checklists from checklist_PA/buffer/out
extent_PA<-extent_multipoly(PA$geometry[[1]]) # lon-min,lon-max,lat-min,lat-max
extent_buffer<-extent_multipoly(buffer$geometry[[1]])
extent_out<-extent_multipoly(out$geometry[[1]])#50km
#chlist=as.data.frame(checklist)
chlist_PA<-subset(checklist_PA,(checklist_PA$lon>=extent_PA$xmin)&(checklist_PA$lon<=extent_PA$xmax)&
                    (checklist_PA$lat>=extent_PA$ymin)&(checklist_PA$lat<=extent_PA$ymax))
chlist_buffer<-subset(checklist_buffer,(checklist_buffer$lon>=extent_buffer$xmin)&(checklist_buffer$lon<=extent_buffer$xmax)&
                        (checklist_buffer$lat>=extent_buffer$ymin)&(checklist_buffer$lat<=extent_buffer$ymax))
chlist_out<-subset(checklist_out,(checklist_out$lon>=extent_out$xmin)&(checklist_out$lon<=extent_out$xmax)&
                     (checklist_out$lat>=extent_out$ymin)&(checklist_out$lat<=extent_out$ymax))

## 4. delete checklists outside the three polygons
###PA
chlist_PA[ , 'inside'] = NA;
for(i in 1:nrow(chlist_PA)){
  chlist_PA$inside[i]<-point.in.multipolygon(chlist_PA$lat[i],chlist_PA$lon[i],PA)
}
chlist_PA<-subset(chlist_PA,chlist_PA$inside!=0)
###buffer
if(nrow(chlist_buffer)!=0){
  chlist_buffer[ , 'inside'] = NA;
  for(i in 1:nrow(chlist_buffer)){
    chlist_buffer$inside[i]<-point.in.multipolygon(chlist_buffer$lat[i],chlist_buffer$lon[i],buffer)
  }
  chlist_buffer<-subset(chlist_buffer,chlist_buffer$inside!=0)
}

###out
chlist_out[ , 'inside'] = NA;
for(i in 1:nrow(chlist_out)){
  chlist_out$inside[i]<-point.in.multipolygon(chlist_out$lat[i],chlist_out$lon[i],out)
}
chlist_out<-subset(chlist_out,chlist_out$inside!=0)

## 5. matching, remove the checklists after matching, save ID
Matching_result<-data.frame(matrix(NA,min(nrow(chlist_PA),nrow(chlist_out)),6))
colnames(Matching_result)<-c("WDPAID","checklist_PA","checklist_buffer","checklist_out","lat","lon")

# 5.1 calculate the mean values of variables
mean_lat<-mean(chlist_PA$z_lat);
mean_lon<-mean(chlist_PA$z_lon);
mean_NPP<-mean(chlist_PA$z_NPP);
mean_Elevation<-mean(chlist_PA$z_Elevation);
mean_TPI<-mean(chlist_PA$z_TPI);
mean_DisToWater<-mean(chlist_PA$z_Distowater);
V_mean<-c(mean_lat, mean_lon, mean_NPP, mean_Elevation, mean_TPI, mean_DisToWater);


## 5.2 checklist_PA>=/<checklist_buffer: if/else
##----- chlist_PA<=chlist_buffer----- matching according to chlist_PA
if(nrow(chlist_PA)<=nrow(chlist_buffer)){
  # 5.3 arrange according to distance to mean values
  chlist_PA[ , 'DisToMean'] = NA;
  for(i in 1:nrow(chlist_PA)){
    V_ch<-c(chlist_PA$z_lat[i], chlist_PA$z_lon[i], chlist_PA$z_NPP[i], chlist_PA$z_Elevation[i], chlist_PA$z_TPI[i], chlist_PA$z_Distowater[i])
    a<-rbind(V_ch,V_mean)
    chlist_PA$DisToMean[i]<-dist(a);
  }
  # arrange according to dist (from highest to lowest)
  chlist_PA <- chlist_PA %>%arrange(desc(DisToMean))
  
  # add dist field to chlist_buffer and chlist_out
  chlist_buffer[ , 'DisToPA'] = NA;
  chlist_out[ , 'DisToPA'] = NA;
  
  # 5.4 calculate distance from the matching checklist
  for(j in 1:nrow(chlist_PA)){
    if(nrow(chlist_out)!=0){
      #chlists outside < PA
      # for each chlist in PA, calculate the dist of buffer from this point
      V_PA<-c(chlist_PA$z_lat[j], chlist_PA$z_lon[j], chlist_PA$z_NPP[j], chlist_PA$z_Elevation[j], chlist_PA$z_TPI[j], chlist_PA$z_Distowater[j])
      for(m in 1:nrow(chlist_buffer)){
        V_buffer<-c(chlist_buffer$z_lat[m],chlist_buffer$z_lon[m],chlist_buffer$z_NPP[m],chlist_buffer$z_Elevation[m],chlist_buffer$z_TPI[m],chlist_buffer$z_Distowater[m])
        b<-rbind(V_PA,V_buffer)
        chlist_buffer$DisToPA[m]<-dist(b);
      }
      for(n in 1:nrow(chlist_out)){
        V_out<-c(chlist_out$z_lat[n],chlist_out$z_lon[n],chlist_out$z_NPP[n],chlist_out$z_Elevation[n],chlist_out$z_TPI[n],chlist_out$z_Distowater[n])
        c<-rbind(V_PA,V_out)
        chlist_out$DisToPA[n]<-dist(c)
      }
      
      # 5.5 arrange checklists according to distance, desc(lat,lon)
      chlist_buffer<-chlist_buffer%>%arrange(DisToPA,desc(z_lat),desc(z_lon))#检查是否先按dis排序，再按lat和lon
      chlist_out<-chlist_out%>%arrange(DisToPA,desc(z_lat),desc(z_lon))
      
      # 5.6 extract the first one, save to results, delete from subsets
      # save
      Matching_result$WDPAID[j]<-PA$WDPAID[1]
      Matching_result$checklist_PA[j]<-chlist_PA$checklist[j]
      Matching_result$checklist_buffer[j]<-chlist_buffer$checklist[1]
      Matching_result$checklist_out[j]<-chlist_out$checklist[1]
      Matching_result$lon[j]<-chlist_PA$lon[j]
      Matching_result$lat[j]<-chlist_PA$lat[j]
      
      #delete
      chlist_buffer<-chlist_buffer[-1, ]
      chlist_out<-chlist_out[-1, ]
    }
  }
}else{
  ##################################################
  ####### chlists in PA > chlists in buffer ########
  ##################################################
  
  ## 5.2 checklist_PA>=/<checklist_buffer: else
  
  #calculate distance of chlist_PA to mean values
  chlist_PA[ , 'DisToMean'] = NA;
  for(i in 1:nrow(chlist_PA)){
    V_ch<-c(chlist_PA$z_lat[i], chlist_PA$z_lon[i], chlist_PA$z_NPP[i], chlist_PA$z_Elevation[i], chlist_PA$z_TPI[i], chlist_PA$z_Distowater[i])
    a<-rbind(V_ch,V_mean)
    chlist_PA$DisToMean[i]<-dist(a);
  }
  
  # 5.3 arrange according to dist of chlist_buffer to mean_buffer(from highest to lowest)
  mean_lat_b<-mean(chlist_buffer$z_lat);
  mean_lon_b<-mean(chlist_buffer$z_lon);
  mean_NPP_b<-mean(chlist_buffer$z_NPP);
  mean_Elevation_b<-mean(chlist_buffer$z_Elevation);
  mean_TPI_b<-mean(chlist_buffer$z_TPI);
  mean_DisToWater_b<-mean(chlist_buffer$z_Distowater);
  V_mean_b<-c(mean_lat_b, mean_lon_b, mean_NPP_b, mean_Elevation_b, mean_TPI_b, mean_DisToWater_b)
  
  if(nrow(chlist_buffer)!=0){
    chlist_buffer[ , 'DisToMean'] = NA;
    for(i in 1:nrow(chlist_buffer)){
      V_buffer<-c(chlist_buffer$z_lat[i], chlist_buffer$z_lon[i], chlist_buffer$z_NPP[i], 
                  chlist_buffer$z_Elevation[i], chlist_buffer$z_TPI[i], chlist_buffer$z_Distowater[i])
      d<-rbind(V_buffer,V_mean_b)#distance to mean in buffer
      chlist_buffer$DisToMean[i]<-dist(d);
    }
    chlist_buffer <- chlist_buffer %>%arrange(desc(DisToMean))
  }
  
  # add dist field to chlist_PA and chlist_out
  chlist_out[ , 'DisToPA'] = NA;
  
  # 5.4 calculate distance from the matching checklist
  if(nrow(chlist_buffer)!=0){
    chlist_PA[ , 'DisTobuffer'] = NA;
    for(j in 1:nrow(chlist_buffer)){
      if(nrow(chlist_out)!=0){
        # for each chlist in buffer, calculate the dist of PA from this point
        V_buffer<-c(chlist_buffer$z_lat[j], chlist_buffer$z_lon[j], chlist_buffer$z_NPP[j], 
                    chlist_buffer$z_Elevation[j], chlist_buffer$z_TPI[j], chlist_buffer$z_Distowater[j])
        for(m in 1:nrow(chlist_PA)){
          V_PA<-c(chlist_PA$z_lat[m],chlist_PA$z_lon[m],chlist_PA$z_NPP[m],
                  chlist_PA$z_Elevation[m],chlist_PA$z_TPI[m],chlist_PA$z_Distowater[m])
          e<-rbind(V_buffer,V_PA)
          chlist_PA$DisTobuffer[m]<-dist(e);
        }
        
        # 5.5 arrange checklists according to distance, desc(lat,lon)
        chlist_PA<-chlist_PA%>%arrange(DisTobuffer,desc(z_lat),desc(z_lon))#先按dis升序排序，再按lat和lon降序
        Matching_result$WDPAID[j]<-PA$WDPAID[1]
        Matching_result$checklist_PA[j]<-chlist_PA$checklist[1]
        Matching_result$checklist_buffer[j]<-chlist_buffer$checklist[j]
        Matching_result$lon[j]<-chlist_PA$lon[1]
        Matching_result$lat[j]<-chlist_PA$lat[1]
        
        # 5.6 extract V_PA of the matching chlist, match out with PA
        V_PA2<-c(chlist_PA$z_lat[1],chlist_PA$z_lon[1],chlist_PA$z_NPP[1],
                 chlist_PA$z_Elevation[1],chlist_PA$z_TPI[1],chlist_PA$z_Distowater[1])
        for(n in 1:nrow(chlist_out)){
          V_out<-c(chlist_out$z_lat[n],chlist_out$z_lon[n],chlist_out$z_NPP[n],
                   chlist_out$z_Elevation[n],chlist_out$z_TPI[n],chlist_out$z_Distowater[n])
          f<-rbind(V_PA2,V_out)
          chlist_out$DisToPA[n]<-dist(f)
        }
        chlist_out<-chlist_out%>%arrange(DisToPA,desc(z_lat),desc(z_lon))
        Matching_result$checklist_out[j]<-chlist_out$checklist[1]
        
        # 5.7 (else) delete the chlist of PA and out
        chlist_PA<-chlist_PA[-1, ]
        chlist_out<-chlist_out[-1, ]
      }
    }}#finish the matching of chlist_buffer
  
  # 5.8 (esle) matching the remaining chlist_PA with out  
  chlist_PA<-chlist_PA%>%arrange(desc(DisToMean))
  for(p in 1:nrow(chlist_PA)){
    if(nrow(chlist_out)!=0){
      V_PA3<-c(chlist_PA$z_lat[p],chlist_PA$z_lon[p],chlist_PA$z_NPP[p],
               chlist_PA$z_Elevation[p],chlist_PA$z_TPI[p],chlist_PA$z_Distowater[p])
      for(q in 1:nrow(chlist_out)){
        V_out2<-c(chlist_out$z_lat[q],chlist_out$z_lon[q],chlist_out$z_NPP[q],
                  chlist_out$z_Elevation[q],chlist_out$z_TPI[q],chlist_out$z_Distowater[q])
        g<-rbind(V_PA3,V_out2)
        chlist_out$DisToPA[q]<-dist(g)
      }
      
      # arrange
      chlist_out<-chlist_out%>%arrange(DisToPA,desc(z_lat),desc(z_lon))
      Matching_result$WDPAID[nrow(chlist_buffer)+p]<-PA$WDPAID[1]
      Matching_result$checklist_PA[nrow(chlist_buffer)+p]<-chlist_PA$checklist[p]
      Matching_result$checklist_out[nrow(chlist_buffer)+p]<-chlist_out$checklist[1]
      Matching_result$lon[nrow(chlist_buffer)+p]<-chlist_PA$lon[p]
      Matching_result$lat[nrow(chlist_buffer)+p]<-chlist_PA$lat[p]
      
      # 5.7 (else) delete the chlist of out
      chlist_out<-chlist_out[-1, ]
    }
  }#finish matching of remaining chlist_PA
  
}#finish else

######################################################################
## all the matching of WDPAID=WDPAID[array_index,1] have finished  ###
########  save all the matching tables of checklist number  ##########
######################################################################

write.csv(Matching_result,paste0(folder_save,"/Matching_result_US",array_index,".csv"),row.names=FALSE)

print("F I N I S H E D")