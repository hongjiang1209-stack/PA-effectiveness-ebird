# 7.Richness_duration_estimation.R
# Purpose:
# Fit nonlinear species richness-duration relationships for matched PA, spillover,
# and control-area checklists for each protected area, and estimate species richness
# under standardized 3-hour observation durations.
#
# Input:
# Combined checklist matching results, filtered checklist datasets for PA, spillover,
# and control-area groups, and a WDPAID list.
#
# Output:
# RDS files containing fitted model parameters, average richness, and estimated
# species richness at 3 hours for each protected area and group.
#
# 
# Notes: Script run on the University of York's HPC, Viking2.



# config ------------------------------------------------------------------

# Set working directory
directory <- '.' 

# Load packages
print("Loading libraries")
###libraries
library(stats)
library(minpack.lm)
library(propagate) #use predictNLS
print('libraries loaded')


time<-"2010"
# time<-"2020"
cargs <- as.numeric(commandArgs(trailingOnly = TRUE))
array_index <- cargs[1]

## predict 1h, 2h, 3h

folder_save <- file.path(directory,paste0("result_",time))

if ( ! dir.exists(folder_save)){
  dir.create(folder_save, showWarnings = FALSE)
}


###############################################
#############  Define a function  #############
###############################################

#object is the nls result (fm)
#x_predict is the point for prediction (duration=60 min in this study)
#no adjustment; interval="confidence"
#refered to 

######## read input datasets ########
##  Matching_results_combined.rds  ##
##  US2020_records_filtered7.rds   ##
##  WDPAID_nofchlist.rds           ##
#####################################
WDPAID<-readRDS(paste0("WDPAID_nofchlist_",time,".rds"))
#WDPAID<-subset(WDPA,(WDPA$nofchlist_PA!=0))
checklist_matching<-readRDS(paste0("Matching_results_combined_",time,".rds"))
checklist_matching<-checklist_matching[,-(5:6)]
checklist_PA<-readRDS(paste0("checklist_PA",time,".rds"))  ############################
checklist_buffer<-readRDS(paste0("checklist_buffer",time,".rds"))   #######################
checklist_out<-readRDS(paste0("checklist_out",time,".rds"))


###########################################
##       Scenario 2: PA+buffer+out       ##
###########################################
checklist_matching<-subset(checklist_matching,!is.na(checklist_matching$checklist_buffer))
group<-c("PA","buffer","out")
result_propagate<-data.frame(matrix(NA,3,10))
colnames(result_propagate)<-c("WDPAID","group","nofchlist","A","B","R_avg","Estim_3H","lwr","upr","Bstart")

j<-array_index
for(i in 1:3){
  line<-i
  ID<-WDPAID[j,1]
  result_propagate$WDPAID[line]<-ID
  result_propagate$group[line]<-group[i]
  matching<-subset(checklist_matching,checklist_matching$WDPAID==ID)
  #matching<-subset(matching,!is.na(matching$checklist_buffer))
  chlist<-matching[,i+1]
  if(i==1){
    sub_chlist<-subset(checklist_PA,checklist_PA$checklist %in% chlist)
  }else if(i==2){
    #chlist<-subset(chlist,!is.na(chlist))
    sub_chlist<-subset(checklist_buffer,checklist_buffer$checklist %in% chlist)
  }else{
    sub_chlist<-subset(checklist_out,checklist_out$checklist %in% chlist)
  }
  
  result_propagate$nofchlist[line]<-nrow(sub_chlist)
  result_propagate$R_avg[line]<-mean(sub_chlist$rich)
  # change duration to hours
  sub_chlist$duration<-sub_chlist$duration/60
  
  if(nrow(sub_chlist)>1){
    rich<-sub_chlist$rich
    duration<-sub_chlist$duration

    re1=tryCatch({
      fm <- nlsLM(rich ~ A*(1-exp(B*duration)),start=list(A=0.9*max(sub_chlist$rich),B=-1))
    },error=function(e){
      "1"
    })
    
    if(length(re1)!=1){
      result_propagate$Bstart[line]<--1
    }else{
      re1=tryCatch({
        fm <- nlsLM(rich ~ A*(1-exp(B*duration)),start=list(A=0.9*max(sub_chlist$rich),B=-0.1))
      },error=function(e){
        "1"
      })
      
      if(length(re1)!=1){
        result_propagate$Bstart[line]<--0.1
      }else{
        re1=tryCatch({
          fm <- nlsLM(rich ~ A*(1-exp(B*duration)),start=list(A=0.9*max(sub_chlist$rich),B=-0.01))
        },error=function(e){
          "1"
        })
        if(length(re1)!=1){
          result_propagate$Bstart[line]<--0.01
        }else{
          re1=tryCatch({
            fm <- nlsLM(rich ~ A*(1-exp(B*duration)),start=list(A=0.9*max(sub_chlist$rich),B=-0.001))
          },error=function(e){
            "1"
          })
          if(length(re1)!=1){
            result_propagate$Bstart[line]<--0.001
          }else{
            re1=tryCatch({
              fm <- nlsLM(rich ~ A*(1-exp(B*duration)),start=list(A=0.9*max(sub_chlist$rich),B=-0.0001))
            },error=function(e){
              "1"
            })
            if(length(re1)!=1){
              result_propagate$Bstart[line]<--0.0001
            }else{
              re1=tryCatch({
                fm <- nlsLM(rich ~ A*(1-exp(B*duration)),start=list(A=0.9*max(sub_chlist$rich),B=-0.0001))
              },error=function(e){
                "1"
              })
              if(length(re1)!=1){
                result_propagate$Bstart[line]<--0.0001
              }else{
                re1=tryCatch({
                  fm <- nlsLM(rich ~ A*(1-exp(B*duration)),start=list(A=0.9*max(sub_chlist$rich),B=0))
                },error=function(e){
                  "1"
                })
                if(length(re1)!=1){
                  result_propagate$Bstart[line]<-0
                }
              }
            }
          }
        }
      }
    }
    
    if(length(re1)!=1){
      result_propagate$A[line]<-summary(fm)[["parameters"]][1,1]
      result_propagate$B[line]<-summary(fm)[["parameters"]][2,1]
      #####
      formula<-stats::formula(fm)
      
      newdata<-data.frame(duration=3)
      result_propagate$Estim_3H[line]<-fm$m$predict(newdata)[1]
      #####
      
      re2=tryCatch({
        z<-predictNLS(fm, data.frame(duration=3), interval = "confidence",alpha = 0.05)
      },error=function(e){
        "1"
      })
      if(length(re2)!=1){
        result_propagate$lwr[line]<-z$summary$`Prop.2.5%`
        result_propagate$upr[line]<-z$summary$`Prop.97.5%`
      }
    }
  }
}
print(paste0(j," in ",time," finished"))


saveRDS(result_propagate,paste0(folder_save,"/richness_duration",time,"_PBO_noW_",array_index,".rds"))
print("F I N I S H E D")
