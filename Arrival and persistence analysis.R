##Arrival & persistence analysis for priority effects review
##Reena Debray

##import libraries
library(ggplot2)
library(readr)
library(RColorBrewer)
library(gridExtra)
getPalette=colorRampPalette(brewer.pal(9,"Set1"))

##import data files: OTU tables ("otu"), taxonomic annotations ("SILVA taxonomy"), trait annotations ("predicted_trait_data")
otus <- read_csv("~/Desktop/Microbiomes as food webs/microbiomepriorityeffects-master/otus.csv")
SILVA_taxonomy <- read_csv("~/Desktop/Microbiomes as food webs/microbiomepriorityeffects-master/SILVA_taxonomy.csv")
predicted_trait_data <- read_csv("~/Desktop/Microbiomes as food webs/microbiomepriorityeffects-master/predicted_trait_data.csv")


##Make a list of OTUs that have predicted traits & are seen in 5 or more subjects
OTU_list<-c()
for (OTU in predicted_trait_data$otu){
  OTU_all_samples<-otus[,c("subject",OTU)]
  
  colnames(OTU_all_samples)<-c("subject","abundance")
  OTU_all_indiv<-OTU_all_samples[OTU_all_samples$abundance>0,]
  
  if (length(unique(OTU_all_indiv$subject))>=5){OTU_list<-c(OTU_list,OTU)}
  
}

##Persistence analysis
summary_all_OTUs<-data.frame(matrix(nrow=0,ncol=26)) #initialize a data frame for mean arrival times, mean persistence, and the arrival-persistence correlation for all the OTUs

for (OTU in OTU_list){
  
  data_per_OTU<-data.frame(matrix(nrow=0,ncol=4)) #initialize a new data frame for each OTU, to calculate the arrival time and persistence in each infant
  
  for (subjectID in levels(factor(otus$subject))){
    
    #subset the data to abundances of one OTU in one subject over time
    otus_subset<-cbind(otus[otus$subject==subjectID,1:4],otus[otus$subject==subjectID,OTU])
    
    #order by time point
    otus_subset$t<-as.numeric(otus_subset$t)
    otus_subset<-otus_subset[order(otus_subset$t),]
    colnames(otus_subset)=c("X1","sampleID","subject","t","relative_abundance")
    
    #calculate arrival time & persistence
    if (!identical(otus_subset[,"relative_abundance"],rep(0,nrow(otus_subset)))){ #if this OTU never appeared in this subject, skip it
      row=1
      while (otus_subset[row,"relative_abundance"]==0){
        row=row+1
      }
      arrival_time<-otus_subset[row,"t"]  #arrival time is the first sample with a non-zero abundance
      post_arrival<-otus_subset[otus_subset$t>arrival_time,"relative_abundance"] #subset the data to the time points after arrival
      persistence<-length(post_arrival[post_arrival!=0])/length(post_arrival) #proportion of samples in which this OTU occurs after arrival
      
      data_per_OTU<-rbind(data_per_OTU,data.frame(OTU,subjectID,arrival_time,persistence))
      
    }
  }
  
  #calculate summary statistics for each OTU across infants
  mean_arrival=mean(data_per_OTU$arrival_time)
  mean_persistence=mean(data_per_OTU$persistence,na.rm=T)
  arrival_persistence_cor=cor.test(data_per_OTU$arrival_time,data_per_OTU$persistence,met="p")$estimate
  arrival_persistence_pvalue=cor.test(data_per_OTU$arrival_time,data_per_OTU$persistence,met="p")$p.value
  num_infants=nrow(data_per_OTU[!is.na(data_per_OTU$persistence),])
  
  #add taxonomy and trait data
  taxa<-SILVA_taxonomy[SILVA_taxonomy$otu==OTU,c(2,4:9)]
  traits<-predicted_trait_data[predicted_trait_data$otu==OTU,2:14]
  
    summary_all_OTUs<-rbind(summary_all_OTUs,data.frame(OTU,mean_arrival,mean_persistence,arrival_persistence_cor,arrival_persistence_pvalue,num_infants,taxa,traits))
}

##summary statistics
#identify early arrivers (1st quartile), late arrivers (4th quartile)
arrival_first_quartile<-summary(summary_all_OTUs$mean_arrival)[2]
arrival_fourth_quartile<-summary(summary_all_OTUs$mean_arrival)[5]
summary_all_OTUs[summary_all_OTUs$mean_arrival<arrival_first_quartile,"arrival_group"]="early"
summary_all_OTUs[summary_all_OTUs$mean_arrival>=arrival_first_quartile & summary_all_OTUs$mean_arrival<=arrival_fourth_quartile,"arrival_group"]="middle"
summary_all_OTUs[summary_all_OTUs$mean_arrival>arrival_fourth_quartile,"arrival_group"]="late"
summary_all_OTUs$arrival_group<-factor(summary_all_OTUs$arrival_group,levels=c("early","middle","late"))

#identify poor persisters (1st quartile), good persisters (4th quartile)
persistence_first_quartile<-summary(summary_all_OTUs$mean_persistence)[2]
persistence_fourth_quartile<-summary(summary_all_OTUs$mean_persistence)[5]
summary_all_OTUs[summary_all_OTUs$mean_persistence<persistence_first_quartile,"persistence_group"]="poor"
summary_all_OTUs[summary_all_OTUs$mean_persistence>=persistence_first_quartile & summary_all_OTUs$mean_persistence<=persistence_fourth_quartile,"persistence_group"]="middle"
summary_all_OTUs[summary_all_OTUs$mean_persistence>persistence_fourth_quartile,"persistence_group"]="good"
summary_all_OTUs$persistence_group<-factor(summary_all_OTUs$persistence_group,levels=c("poor","middle","good"))

#identify time-dependent colonizers
summary_all_OTUs[!is.na(summary_all_OTUs$arrival_persistence_pvalue) & summary_all_OTUs$arrival_persistence_pvalue<0.05 & summary_all_OTUs$arrival_persistence_cor>0,"priority_effects"]="Prefers late arrival"
summary_all_OTUs[!is.na(summary_all_OTUs$arrival_persistence_pvalue) & summary_all_OTUs$arrival_persistence_pvalue<0.05 & summary_all_OTUs$arrival_persistence_cor<0,"priority_effects"]="Prefers early arrival"
summary_all_OTUs[!is.na(summary_all_OTUs$arrival_persistence_pvalue) & summary_all_OTUs$arrival_persistence_pvalue>=0.05,"priority_effects"]="None"

#plot time-dependent colonizers by phylum and class
ggplot(summary_all_OTUs[summary_all_OTUs$priority_effects!="None" & !is.na(summary_all_OTUs$priority_effects),],aes(priority_effects,fill=Phylum))+geom_bar()+theme_bw()+xlab("")+scale_fill_brewer(palette="Set1")
colourCount = length(unique(summary_all_OTUs$Class))
ggplot(summary_all_OTUs[summary_all_OTUs$priority_effects!="None" & !is.na(summary_all_OTUs$priority_effects),],aes(priority_effects,fill=Class))+geom_bar()+theme_bw()+xlab("")+scale_fill_manual(values = getPalette(colourCount))


##Pairwise analysis
#start with the OTUs that do better if they arrive early, and determine all possible pairs
#NOTE that this list does not take multiple testing into account
prefer_early<-as.character(summary_all_OTUs[!is.na(summary_all_OTUs$priority_effects) & summary_all_OTUs$priority_effects=="Prefers early arrival","OTU"]) #names of OTUs in the prefer-early-arrival group
prefer_early_pairs<-combn(prefer_early,2) #all possible pairs of early-preferring OTUs

#for each pair of OTUs, test whether the relative order of arrival is a BETTER predictor than chronological time.

summary_all_pairs<-data.frame(matrix(nrow=0,ncol=17))
for (i in seq(1,ncol(prefer_early_pairs))){
  OTU_A<-prefer_early_pairs[,i][1]
  OTU_B<-prefer_early_pairs[,i][2]
  
  data_per_pair<-data.frame(matrix(nrow=0,ncol=6))
  
  #loop through subjects
  #in each subject, determine which OTU arrived first
  for (subjectID in levels(factor(otus$subject))){
    
    #subset the data to abundances of two focal OTUS in one subject over time
    otus_subset<-cbind(otus[otus$subject==subjectID,1:4],otus[otus$subject==subjectID,OTU_A],otus[otus$subject==subjectID,OTU_B])
    otus_subset$t<-as.numeric(otus_subset$t)
    otus_subset<-otus_subset[order(otus_subset$t),]
    colnames(otus_subset)=c("X1","sampleID","subject","t","OTU_A","OTU_B")
    
    #check that both OTUs appear at least once before the last sample -- otherwise, skip this subject
    if (!identical(otus_subset[-c(nrow(otus_subset)),"OTU_A"],rep(0,nrow(otus_subset)-1)) & !identical(otus_subset[-c(nrow(otus_subset)),"OTU_B"],rep(0,nrow(otus_subset)-1))){
      
      #calculate arrival time & persistence of A
      row=1
      while (otus_subset[row,"OTU_A"]==0){row=row+1}
      arrival_time_A<-otus_subset[row,"t"]
      post_arrival_A<-otus_subset[otus_subset$t>arrival_time_A,"OTU_A"] #relative abundance of A after it arrives
      persistence_A<-length(post_arrival_A[post_arrival_A!=0])/length(post_arrival_A) #proportion of samples in which this OTU occurs after arrival
      
      #calculate arrival time & persistence of B
      row=1
      while (otus_subset[row,"OTU_B"]==0){row=row+1}
      arrival_time_B<-otus_subset[row,"t"]
      post_arrival_B<-otus_subset[otus_subset$t>arrival_time_B,"OTU_B"] #relative abundance of B after both arrive
      persistence_B<-length(post_arrival_B[post_arrival_B!=0])/length(post_arrival_B) #proportion of samples in which this OTU occurs after arrival

      arrival_difference<-arrival_time_A-arrival_time_B #if A arrives first, this is negative
      co_occurrence<-nrow(otus_subset[otus_subset$OTU_A!=0 & otus_subset$OTU_B!=0,])/nrow(otus_subset[otus_subset$OTU_A!=0 | otus_subset$OTU_B!=0,]) #This is a measure of the average overlap between OTUs within infants. The value is lowest if they tend to colonize the same infants but do not coexist at the same time.
    
      
      data_per_pair<-rbind(data_per_pair,data.frame(arrival_time_A,arrival_time_B,arrival_difference,persistence_A,persistence_B,co_occurrence))
    }
  }
  #check whether we have a roughly similar distribution of arrival times between the two OTUs. If one always arrives first, then we can't look for priority effects.
  if (nrow(data_per_pair)>=5 & nrow(data_per_pair[data_per_pair$arrival_difference<=0,])/nrow(data_per_pair)<0.75 & nrow(data_per_pair[data_per_pair$arrival_difference>=0,])/nrow(data_per_pair)<0.75){
      chrono_R2_A<-summary(lm(persistence_A~arrival_time_A,data_per_pair))$r.squared #how well does chronological time explain the persistence of A?
      order_R2_A<-summary(lm(persistence_A~arrival_difference,data_per_pair))$r.squared  #how well does arrival order explain the persistence of A?
      order_directionality_A<-summary(lm(persistence_A~arrival_difference,data_per_pair))$coef[2,1] #if there is a NEGATIVE relationship, that means A persists better when it arrives before B
        
        
      chrono_R2_B<-summary(lm(persistence_B~arrival_time_B,data_per_pair))$r.squared #how well does chronological time explain the persistence of B?
      order_R2_B<-summary(lm(persistence_B~arrival_difference,data_per_pair))$r.squared  #how well does arrival order explain the persistence of B?
      order_directionality_B<-summary(lm(persistence_B~arrival_difference,data_per_pair))$coef[2,1] #if there is a POSITIVE relationship, that means B persists better when it arrives before A
     
      avg_co_occurrence<-mean(data_per_pair$co_occurrence) 
      num_infants<-nrow(data_per_pair) #total number of infant samples where both are observed at least once
      
      SILVA_A<-summary_all_OTUs[summary_all_OTUs$OTU==OTU_A,8:11]
      SILVA_B<-summary_all_OTUs[summary_all_OTUs$OTU==OTU_B,8:11]

      
      summary_all_pairs<-rbind(summary_all_pairs,data.frame(OTU_A,OTU_B,chrono_R2_A,order_R2_A,order_directionality_A,chrono_R2_B,order_R2_B,order_directionality_B,avg_co_occurrence,num_infants,SILVA_A,SILVA_B))
  }
}

#using taxonomy information, determine whether the two OTUs are from the same or different groups
summary_all_pairs$same_phylum<-summary_all_pairs$Phylum==summary_all_pairs$Phylum.1
summary_all_pairs$same_class<-summary_all_pairs$Class==summary_all_pairs$Class.1
summary_all_pairs$same_order<-summary_all_pairs$Order==summary_all_pairs$Order.1
summary_all_pairs$same_family<-summary_all_pairs$Family==summary_all_pairs$Family.1

#bin OTU pairs according to whether chronological time is a better predictor of persistence than order for both of them ("Arrival time"), arrival order is better than chrono time for one of the two ("A or B"), and arrival order is better than chrono time for both OTUs ("Both A and B")
summary_all_pairs<-summary_all_pairs[!is.na(summary_all_pairs$order_R2_A & summary_all_pairs$order_R2_B),]

summary_all_pairs[summary_all_pairs$order_R2_A>summary_all_pairs$chrono_R2_A & summary_all_pairs$order_directionality_A<0 & summary_all_pairs$order_R2_B>summary_all_pairs$chrono_R2_B & summary_all_pairs$order_directionality_B>0,"preemption"]="Both A and B"

summary_all_pairs[(summary_all_pairs$order_R2_A>summary_all_pairs$chrono_R2_A & summary_all_pairs$order_directionality_A<0) & (summary_all_pairs$order_R2_B<=summary_all_pairs$chrono_R2_B | summary_all_pairs$order_directionality_B<0),"preemption"]="A or B"

summary_all_pairs[(summary_all_pairs$order_R2_A<=summary_all_pairs$chrono_R2_A | summary_all_pairs$order_directionality_A>0) & (summary_all_pairs$order_R2_B>summary_all_pairs$chrono_R2_B & summary_all_pairs$order_directionality_B>0),"preemption"]="A or B"

summary_all_pairs[(summary_all_pairs$order_R2_A<=summary_all_pairs$chrono_R2_A & summary_all_pairs$order_R2_B<=summary_all_pairs$chrono_R2_B) | (summary_all_pairs$order_directionality_A>0 | summary_all_pairs$order_directionality_B<0),"preemption"]="Arrival time"

summary_all_pairs$preemption<-factor(summary_all_pairs$preemption,levels=c("Arrival time","A or B","Both A and B"))

#plot the strength of these pairwise relationships according to co-occurrence rates
ggplot(summary_all_pairs,aes(avg_co_occurrence,fill=preemption))+geom_density(alpha=0.5)+facet_wrap(~preemption,nrow=3)+theme_bw()+guides(fill=F)+xlab("Co-occurrence (number of infants)")

#plot the strength of these pairwise relationships according to taxonomic relatedness
grid.arrange(ggplot(summary_all_pairs,aes(preemption,fill=same_phylum))+geom_bar(position="dodge")+theme_bw()+scale_fill_brewer(palette="Set2")+xlab(""),ggplot(summary_all_pairs,aes(preemption,fill=same_class))+geom_bar(position="dodge")+theme_bw()+scale_fill_brewer(palette="Set2")+xlab(""))



