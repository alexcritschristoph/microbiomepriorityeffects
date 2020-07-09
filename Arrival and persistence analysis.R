##Arrival & persistence analysis for priority effects review
##Reena Debray

##import libraries
library(ggplot2)
library(readr)
library(RColorBrewer)
library(gridExtra)
library(ape)
library(vegan)
getPalette=colorRampPalette(brewer.pal(9,"Set1"))

##import OTU tables and OTU taxonomy
#human data: OTU tables ("humangut_otus"), taxonomic annotations ("humangut_taxonomy"), trait annotations ("predicted_trait_data")
humangut_otus <- read_csv("~/Desktop/Microbiomes as food webs/Human gut data/otus.csv")
humangut_taxonomy <- read_csv("~/Desktop/Microbiomes as food webs/Human gut data/SILVA_taxonomy.csv")
predicted_trait_data <- read_csv("~/Desktop/Microbiomes as food webs/Human gut data/predicted_trait_data.csv")
#mouse data: OTU tables ("murine_otus"), taxonomic annotations ("murine_taxonomy")
murine_otus<-read_csv("~/Desktop/Microbiomes as food webs/Mouse gut data/murine_OTUs.csv")
murine_taxonomy <- read_excel("~/Desktop/Microbiomes as food webs/Mouse gut data/murine_taxonomy.xlsx")
#cow data: OTU tables ("rumen_otus"), taxonomic annotations ("rumen_taxonomy")
rumen_otus <- read_csv("~/Desktop/Microbiomes as food webs/Rumen gut data/rumen_OTUs.csv")
rumen_taxonomy <- read_csv("~/Desktop/Microbiomes as food webs/Rumen gut data/Rumen gut taxonomy.csv")



------------------------------------------------------------------------------------------------------



##This function creates a data frame of the relative abundances of OTU X in subject Y over time
df_subset<-function(OTU_table,OTU_name,subject_name){
  otus_subset<-OTU_table[OTU_table$subject==subject_name,c("t","subject","sampleID",OTU_name)]
  otus_subset$t<-as.numeric(otus_subset$t)
  otus_subset<-otus_subset[order(otus_subset$t),] #order by time point
 
  return(otus_subset)
}

##This function returns the first time at which relative abundance of OTU X in subject Y is > 0 (or NA if this never occurred)
calculate_arrival_time<-function(otus_subset,OTU_name){
  if (sum(otus_subset[,OTU_name]>0)){ #if this OTU never appeared in this subject, skip it
  i=1
  while (otus_subset[i,OTU_name]==0){i=i+1}
  return(as.numeric(otus_subset[i,"t"])) #arrival time is the first sample with a non-zero abundance
  }
  else {return(NA)}
}

##This function returns persistence (the proportion of a certain window (set by obs_length) after arrival in which OTU X was detected)
calculate_persistence<-function(otus_subset,OTU_name,arrival_time,obs_length){
  
  if (is.na(arrival_time)) {return(NA)}
  if (otus_subset[nrow(otus_subset),"t"]-arrival_time<obs_length) {return(NA)} #if OTU arrived too close to the end of the study, return NA
  if ((otus_subset[match(arrival_time,otus_subset$t)+1,"t"]-arrival_time)>obs_length) {return(NA)} #if no samples were taken within the observation window, return NA
  
  else {
    start_time<-arrival_time
    end_time<-arrival_time+obs_length
    end_i<-match(as.numeric(tail(otus_subset[otus_subset$t<=end_time,"t"],1)),otus_subset$t)
    persistence=0
    i=match(start_time,otus_subset$t)+1
    
    while (i<end_i){
    
        while (otus_subset[i,OTU_name]>0 & i<end_i){i=i+1} #find next dropout (abundance=0) event
        dropout_time<-otus_subset[i-1,"t"] #last time it was observed
        persistence=persistence+(dropout_time-start_time)

        while (otus_subset[i,OTU_name]==0 & i<end_i){i=i+1} #find the next colonization event
        start_time<-otus_subset[i,"t"]
        i=match(start_time,otus_subset$t)
  }
  
  if (otus_subset[end_i,OTU_name]>0){persistence=persistence+otus_subset[end_i,"t"]-otus_subset[end_i-1,"t"]}
    
  obs_window<-otus_subset[end_i,"t"]-(arrival_time) #actual length of time observed (may be slightly different from obs_length depending on when samples were taken)
  prop_persistence<-as.numeric(persistence/obs_window)
    
  return(prop_persistence)
  }
}



------------------------------------------------------------------------------------------------------




##Within each dataset, make a list of OTUs that are seen in at least 20% of hosts
#human gut
humangut_OTU_list<-c()
for (OTU in colnames(humangut_otus)[5:2420]){ #loop through all OTUs detected at any point in the study (2416 OTUs)
  OTU_all_samples<-humangut_otus[,c("subject",OTU)] #subset to this specific OTU
  colnames(OTU_all_samples)<-c("subject","abundance")
  OTU_all_indiv<-OTU_all_samples[OTU_all_samples$abundance>0,] #subset to samples in which this OTU had non-zero abundance
  if (length(unique(OTU_all_indiv$subject))>=12){humangut_OTU_list<-c(humangut_OTU_list,OTU)} #count the number of unique infant hosts in which this OTU had non-zero abundance for any amount of time
  #consider this OTU only if it is present in 20% or more of individuals (11.2)
  #Final list includes 670 OTUs
}
#mouse gut
murine_OTU_list<-c()
for (OTU in colnames(murine_otus[5:3099])){ #loop through all OTUs detected at any point in the study (3095)
  OTU_all_samples<-murine_otus[,c("subject",OTU)] #subset to this specific OTU
  colnames(OTU_all_samples)<-c("subject","abundance")
  OTU_all_indiv<-OTU_all_samples[OTU_all_samples$abundance>0,] #subset to samples in which this OTU had non-zero abundance
  if (length(unique(OTU_all_indiv$subject))>=3){murine_OTU_list<-c(murine_OTU_list,OTU)} #count the number of unique individuals in which this OTU had non-zero abundance for any amount of time
  #consider this OTU only if it is present in 20% or more of individuals (2.4)
  #Final list includes 971 OTUs
}
murine_otus<-murine_otus[murine_otus$t<300,] #remove the last time-point because there was a 6-month gap in sampling
#cow gut
rumen_OTU_list<-c()
for (OTU in colnames(rumen_otus[4:2547])){ #loop through all OTUs detected at any point in the study (2544)
  OTU_all_samples<-rumen_otus[,c("subject",OTU)] #subset to this specific OTU
  colnames(OTU_all_samples)<-c("subject","abundance")
  OTU_all_indiv<-OTU_all_samples[OTU_all_samples$abundance>0,] #subset to samples in which this OTU had non-zero abundance
  if (length(unique(OTU_all_indiv$subject))>=9){rumen_OTU_list<-c(rumen_OTU_list,OTU)} #count the number of unique individuals in which this OTU had non-zero abundance for any amount of time
  #consider this OTU only if it is present in 20% or more of individuals (9)
  #Final list includes 2544 OTUs (because the original study had already subsetted to OTUs present in 80% of individuals)
  }



------------------------------------------------------------------------------------------------------


##identify OTUs that are sensitive to microbiome composition upon arrival (t=0)
#human gut
summary_humangut_OTUs<-data.frame(matrix(nrow=0,ncol=13))
humangut_list<-list()
i=1

for (OTU in humangut_OTU_list){
  data_per_OTU<-data.frame(matrix(nrow=0,ncol=4))
  
  for (subjectID in levels(factor(otus$subject))){
    otus_subset<-df_subset(humangut_otus,OTU,subjectID)
    arrival_time<-calculate_arrival_time(otus_subset,OTU)
    persistence<-calculate_persistence(otus_subset,OTU,arrival_time,6)
    data_per_OTU<-rbind(data_per_OTU,data.frame(OTU,subjectID,arrival_time,persistence))
  }
  
  mean_arrival<-mean(data_per_OTU$arrival_time,na.rm=T)
  mean_persistence<-mean(data_per_OTU$persistence,na.rm=T)
  num_hosts<-nrow(data_per_OTU[!is.na(data_per_OTU$persistence),])
  taxonomy<-humangut_taxonomy[humangut_taxonomy$otu==OTU,c(2,4:9)]
  data_per_OTU<-data_per_OTU[!is.na(data_per_OTU$arrival_time) & !is.na(data_per_OTU$persistence),]
  
   if (num_hosts>=12 & sd(data_per_OTU$persistence)>0){
  
      current_comm_analysis<-data.frame(matrix(nrow=12,ncol=2419))
      colnames(current_comm_analysis)=c("subject","arrival_time","persistence",colnames(otus[5:2420]))
      for (row in seq(1,nrow(data_per_OTU))){
          subject<-as.character(data_per_OTU[row,"subjectID"])
          arrival_time<-data_per_OTU[row,"arrival_time"]
          persistence<-data_per_OTU[row,"persistence"]
  
          #extract current community, remove abundance of focal OTU
          current_community<-data.frame(otus[otus$subject==subject & otus$t==arrival_time,5:2420])
    
          current_comm_analysis[row,1:3]<-c(subject,arrival_time,persistence)
          current_comm_analysis[row,4:2419]<-current_community
          }
        current_comm_analysis$persistence<-as.numeric(current_comm_analysis$persistence)
        focal<-match(OTU,colnames(current_comm_analysis))
        current_comm_analysis<-current_comm_analysis[,-c(focal)]
  
      dist<-vegdist(current_comm_analysis[,4:2418], method="bray")
      a<-adonis(dist~persistence,permutations=1000,data=current_comm_analysis)
      persistence_F<-a$aov.tab[1,4]
      persistence_pvalue<-a$aov.tab[1,6]
  

      summary_humangut_OTUs<-rbind(summary_humangut_OTUs,data.frame(OTU,mean_arrival,mean_persistence,num_hosts,taxonomy,persistence_F,persistence_pvalue))
      humangut_list[[i]]<-data_per_OTU
      i=i+1
     }  
  }
  

#mouse gut
summary_murine_OTUs<-data.frame(matrix(nrow=0,ncol=13)) #initialize a data frame for mean arrival times, mean persistence, and the arrival-persistence correlation for all the OTUs
murine_list<-list()
i=1

for (OTU in murine_OTU_list){
  data_per_OTU<-data.frame(matrix(nrow=0,ncol=4))
  
  for (subjectID in levels(factor(murine_otus$subject))){
    otus_subset<-df_subset(murine_otus,OTU,subjectID)
    arrival_time<-calculate_arrival_time(otus_subset,OTU)
    persistence<-calculate_persistence(otus_subset,OTU,arrival_time,30) #1 month rather than 6 months?
    data_per_OTU<-rbind(data_per_OTU,data.frame(OTU,subjectID,arrival_time,persistence))
  }
  
  mean_arrival<-mean(data_per_OTU$arrival_time,na.rm=T)
  mean_persistence<-mean(data_per_OTU$persistence,na.rm=T)
  num_hosts<-nrow(data_per_OTU[!is.na(data_per_OTU$persistence),])
  taxonomy<-murine_taxonomy[murine_taxonomy$OTU==OTU,c(5,7,9,11,13)]
  data_per_OTU<-data_per_OTU[!is.na(data_per_OTU$arrival_time) & !is.na(data_per_OTU$persistence),]

  if (num_hosts>=3 & sd(data_per_OTU$persistence)>0){
      current_comm_analysis<-data.frame(matrix(nrow=0,ncol=3098))
      colnames(current_comm_analysis)=c("subject","arrival_time","persistence",colnames(murine_otus[5:3099]))
      for (row in seq(1,nrow(data_per_OTU))){
          subject<-as.character(data_per_OTU[row,"subjectID"])
          arrival_time<-data_per_OTU[row,"arrival_time"]
          persistence<-data_per_OTU[row,"persistence"]
  
          #extract current community, remove abundance of focal OTU
          current_community<-data.frame(murine_otus[murine_otus$subject==subject & murine_otus$t==arrival_time,5:3099])
    
          current_comm_analysis[row,1:3]<-c(subject,arrival_time,persistence)
          current_comm_analysis[row,4:3098]<-current_community
          }
        current_comm_analysis$persistence<-as.numeric(current_comm_analysis$persistence)
        focal<-match(OTU,colnames(current_comm_analysis))
        current_comm_analysis<-current_comm_analysis[,-c(focal)]
  
      dist<-vegdist(current_comm_analysis[,3:3097], method="bray")
      a<-adonis(dist~persistence,permutations=1000,data=current_comm_analysis)
      persistence_F<-a$aov.tab[1,4]
      persistence_pvalue<-a$aov.tab[1,6]
  
    
      summary_murine_OTUs<-rbind(summary_murine_OTUs,data.frame(OTU,mean_arrival,mean_persistence,num_hosts,taxonomy,persistence_F,persistence_pvalue))
      murine_list[[i]]<-data_per_OTU
      i=i+1
  }
 
}



#cow gut

summary_rumen_OTUs<-data.frame(matrix(nrow=0,ncol=13)) #initialize a data frame for mean arrival times, mean persistence, and the arrival-persistence correlation for all the OTUs
rumen_list<-list()
i=1

for (OTU in rumen_OTU_list){
  data_per_OTU<-data.frame(matrix(nrow=0,ncol=4))
  
  for (subjectID in levels(factor(rumen_otus$subject))){
    otus_subset<-df_subset(rumen_otus,OTU,subjectID)
    arrival_time<-calculate_arrival_time(otus_subset,OTU)
    persistence<-calculate_persistence(otus_subset,OTU,arrival_time,180)
    data_per_OTU<-rbind(data_per_OTU,data.frame(OTU,subjectID,arrival_time,persistence))
  }
  
  mean_arrival<-mean(data_per_OTU$arrival_time,na.rm=T)
  mean_persistence<-mean(data_per_OTU$persistence,na.rm=T)
  num_hosts<-nrow(data_per_OTU[!is.na(data_per_OTU$persistence),])
  
  taxonomy<-rumen_taxonomy[rumen_taxonomy$X1==OTU,3:8]
  phylum<-substr(taxonomy[1],4,nchar(taxonomy[1]))
  class<-substr(taxonomy[2],4,nchar(taxonomy[2]))
  order<-substr(taxonomy[3],4,nchar(taxonomy[3]))
  family<-substr(taxonomy[4],4,nchar(taxonomy[4]))
  genus<-substr(taxonomy[5],4,nchar(taxonomy[5]))
  species<-substr(taxonomy[6],4,nchar(taxonomy[6]))
  
  data_per_OTU<-data_per_OTU[!is.na(data_per_OTU$arrival_time) & !is.na(data_per_OTU$persistence),]
  
    if (num_hosts>=9 & sd(data_per_OTU$persistence)>0){
      current_comm_analysis<-data.frame(matrix(nrow=0,ncol=2547))
      colnames(current_comm_analysis)=c("subject","arrival_time","persistence",colnames(rumen_otus[4:2547]))
      for (row in seq(1,nrow(data_per_OTU))){
          subject<-as.character(data_per_OTU[row,"subjectID"])
          arrival_time<-data_per_OTU[row,"arrival_time"]
          persistence<-data_per_OTU[row,"persistence"]
  
          #extract current community, remove abundance of focal OTU
          current_community<-data.frame(rumen_otus[rumen_otus$subject==subject & rumen_otus$t==arrival_time,5:2547])
    
          current_comm_analysis[row,1:3]<-c(subject,arrival_time,persistence)
          current_comm_analysis[row,4:2547]<-current_community
          }
        current_comm_analysis$persistence<-as.numeric(current_comm_analysis$persistence)
        focal<-match(OTU,colnames(current_comm_analysis))
        current_comm_analysis<-current_comm_analysis[,-c(focal)]
  
      dist<-vegdist(current_comm_analysis[,3:2546], method="bray")
      a<-adonis(dist~persistence,permutations=1000,data=current_comm_analysis)
      persistence_F<-a$aov.tab[1,4]
      persistence_pvalue<-a$aov.tab[1,6]
  
    
      summary_rumen_OTUs<-rbind(summary_rumen_OTUs,data.frame(OTU,mean_arrival,mean_persistence,num_hosts,taxonomy,persistence_F,persistence_pvalue))
      rumen_list[[i]]<-data_per_OTU
      i=i+1
    }
  }

