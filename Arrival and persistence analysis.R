##Arrival & persistence analysis for priority effects review
##Reena Debray

##import data files: OTU tables ("otu"), taxonomic annotations ("SILVA taxonomy"), trait annotations ("predicted_trait_data")
otus <- read_csv("~/Desktop/Microbiomes as food webs/microbiomepriorityeffects-master/otus.csv")
SILVA_taxonomy <- read_csv("~/Desktop/Microbiomes as food webs/microbiomepriorityeffects-master/SILVA_taxonomy.csv")
predicted_trait_data <- read_csv("~/Desktop/Microbiomes as food webs/microbiomepriorityeffects-master/predicted_trait_data.csv")


#Make a list of OTUs that have predicted traits & are seen in 5 or more subjects
OTU_list<-c()
for (OTU in predicted_trait_data$otu){
  OTU_all_samples<-otus[,c("subject",OTU)]
  
  colnames(OTU_all_samples)<-c("subject","abundance")
  OTU_all_indiv<-OTU_all_samples[OTU_all_samples$abundance>0,]
  
  if (length(unique(OTU_all_indiv$subject))>=5){OTU_list<-c(OTU_list,OTU)}
  
}

#Persistence analysis
indiv_OTU_data<-vector("list",length(OTU_list))
i=1
summary_all_OTUs<-data.frame(matrix(nrow=0,ncol=25))

for (OTU in OTU_list){
  
  data_per_OTU<-data.frame(matrix(nrow=0,ncol=4))
  
  for (subjectID in levels(factor(otus$subject))){
    
    #subset the data to abundances of one OTU in one subject over time
    otus_subset<-cbind(otus[otus$subject==subjectID,1:4],otus[otus$subject==subjectID,OTU])
    
    #order by time point
    otus_subset$t<-as.numeric(otus_subset$t)
    otus_subset<-otus_subset[order(otus_subset$t),]
    colnames(otus_subset)=c("X1","sampleID","subject","t","relative_abundance")
    
    #calculate arrival time & persistence
    if (!identical(otus_subset[,"relative_abundance"],rep(0,nrow(otus_subset)))){
      row=1
      while (otus_subset[row,"relative_abundance"]==0){
        row=row+1
      }
      arrival_time<-otus_subset[row,"t"] #first non-zero abundance
      post_arrival<-otus_subset[otus_subset$t>arrival_time,"relative_abundance"] #relative abundances after arrival
      persistence<-length(post_arrival[post_arrival!=0])/length(post_arrival) #proportion of samples in which this OTU occurs after arrival
      
      data_per_OTU<-rbind(data_per_OTU,data.frame(OTU,subjectID,arrival_time,persistence))
      
    }
  }
  
  #calculate summary statistics for each OTU across infants
  mean_arrival=mean(data_per_OTU$arrival_time)
  mean_persistence=mean(data_per_OTU$persistence,na.rm=T)
  arrival_persistence_cor=cor.test(data_per_OTU$arrival_time,data_per_OTU$persistence,met="s")$estimate
  arrival_persistence_pvalue=cor.test(data_per_OTU$arrival_time,data_per_OTU$persistence,met="s")$p.value
  
  #add taxonomy and trait data
  taxa<-SILVA_taxonomy[SILVA_taxonomy$otu==OTU,c(2,4:9)]
  traits<-predicted_trait_data[predicted_trait_data$otu==OTU,2:14]
  
  summary_all_OTUs<-rbind(summary_all_OTUs,data.frame(OTU,mean_arrival,mean_persistence,arrival_persistence_cor,arrival_persistence_pvalue,taxa,traits))
  indiv_OTU_data[i]=data_per_OTU
  i=i+1
}

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

