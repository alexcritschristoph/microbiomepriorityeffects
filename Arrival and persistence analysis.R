##Arrival & persistence analysis for priority effects review
##Reena Debray

##import libraries
library(ggplot2)
library(readr)
library(RColorBrewer)
library(gridExtra)
getPalette=colorRampPalette(brewer.pal(9,"Set1"))

##import data files: OTU tables ("humangut_otus"), taxonomic annotations ("humangut_taxonomy"), trait annotations ("predicted_trait_data")
humangut_otus <- read_csv("~/Desktop/Microbiomes as food webs/Human gut data/otus.csv")
humangut_taxonomy <- read_csv("~/Desktop/Microbiomes as food webs/Human gut data/SILVA_taxonomy.csv")
predicted_trait_data <- read_csv("~/Desktop/Microbiomes as food webs/Human gut data/predicted_trait_data.csv")

##Make a list of OTUs that are seen in 20% or more of individuals
humangut_OTU_list<-c()
for (OTU in colnames(humangut_otus)[5:2420]){ #loop through all OTUs detected at any point in the study (2416 unique OTUs)
  OTU_all_samples<-humangut_otus[,c("subject",OTU)] #subset to this specific OTU
  colnames(OTU_all_samples)<-c("subject","abundance")
  OTU_all_indiv<-OTU_all_samples[OTU_all_samples$abundance>0,] #subset to samples in which this OTU had non-zero abundance
  if (length(unique(OTU_all_indiv$subject))>=12){humangut_OTU_list<-c(humangut_OTU_list,OTU)} #count the number of unique infant hosts in which this OTU had non-zero abundance for any amount of time
  #consider this OTU only if it is present in 20% or more of individuals (11.2)
}
#We are considering 670 OTUs in this analysis

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

#identify time-dependent colonizers
summary_all_OTUs[!is.na(summary_all_OTUs$arrival_persistence_pvalue) & summary_all_OTUs$arrival_persistence_pvalue<0.05 & summary_all_OTUs$arrival_persistence_cor>0,"priority_effects"]="Prefers late arrival"
summary_all_OTUs[!is.na(summary_all_OTUs$arrival_persistence_pvalue) & summary_all_OTUs$arrival_persistence_pvalue<0.05 & summary_all_OTUs$arrival_persistence_cor<0,"priority_effects"]="Prefers early arrival"
summary_all_OTUs[!is.na(summary_all_OTUs$arrival_persistence_pvalue) & summary_all_OTUs$arrival_persistence_pvalue>=0.05,"priority_effects"]="None"

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

#for pairs in which one OTU preempts the other but not the other way around -- fill in their persistence and occurrence values to look for a difference
one_way_preemption<-summary_all_pairs[summary_all_pairs$preemption=="A or B",]
for (i in seq(1,nrow(one_way_preemption))){
  
  #identify pairs in which A is affected by arrival order but not B
  #(within a pair, whether each OTU is A or B is random)
  if (one_way_preemption[i,"order_R2_A"]>one_way_preemption[i,"chrono_R2_A"] & one_way_preemption[i,"order_R2_B"]<=one_way_preemption[i,"chrono_R2_B"] ){
    focal<-as.character(one_way_preemption[i,"OTU_A"])
    one_way_preemption[i,"focal"]=focal #this is the OTU that IS affected by arrival order
    one_way_preemption[i,"focal_persistence"]=summary_all_OTUs[summary_all_OTUs$OTU==focal,"mean_persistence"]
    one_way_preemption[i,"focal_occurrence"]=summary_all_OTUs[summary_all_OTUs$OTU==focal,"num_infants"]
    
    non_focal<-as.character(one_way_preemption[i,"OTU_B"])
    one_way_preemption[i,"non_focal"]=non_focal #this is the OTU that is NOT affected by arrival order
    one_way_preemption[i,"non_focal_persistence"]=summary_all_OTUs[summary_all_OTUs$OTU==non_focal,"mean_persistence"]
    one_way_preemption[i,"non_focal_occurrence"]=summary_all_OTUs[summary_all_OTUs$OTU==non_focal,"num_infants"]
  }
   #identify pairs in which B is affected by arrival order but not A
    if (one_way_preemption[i,"order_R2_A"]<=one_way_preemption[i,"chrono_R2_A"] & one_way_preemption[i,"order_R2_B"]>one_way_preemption[i,"chrono_R2_B"] ){
      focal<-as.character(one_way_preemption[i,"OTU_B"])
      one_way_preemption[i,"focal"]=focal #this is the OTU that IS affected by arrival order
      one_way_preemption[i,"focal_persistence"]=summary_all_OTUs[summary_all_OTUs$OTU==focal,"mean_persistence"]
      one_way_preemption[i,"focal_occurrence"]=summary_all_OTUs[summary_all_OTUs$OTU==focal,"num_infants"]
    
        non_focal<-as.character(one_way_preemption[i,"OTU_A"])
      one_way_preemption[i,"non_focal"]=non_focal #this is the OTU that is NOT affected by arrival order
      one_way_preemption[i,"non_focal_persistence"]=summary_all_OTUs[summary_all_OTUs$OTU==non_focal,"mean_persistence"]
      one_way_preemption[i,"non_focal_occurrence"]=summary_all_OTUs[summary_all_OTUs$OTU==non_focal,"num_infants"]
    }
}

##Ecological network analysis
#run spiec-easi and extract edge weights
tmp<-as.matrix(data.frame(otus[,OTU_list]))
se.mb.amgut <- spiec.easi(tmp, method='mb', lambda.min.ratio=1e-2,nlambda=20, pulsar.params=list(rep.num=50))
ig.mb     <- adj2igraph(getRefit(se.mb.amgut))
plot(ig.mb,vertex.size=1, vertex.label=NA)
sebeta <- symBeta(getOptBeta(se.mb.amgut), mode='ave')
elist.mb <- Matrix::summary(sebeta)


##Code for figures in document
#figure 1 - arrival time vs persistence across OTUs
ggplot(summary_all_OTUs,aes(mean_arrival,mean_persistence,fill=num_infants))+geom_point(shape=21,size=2)+scale_fill_gradient(low="red",high="yellow",name="Population\nfrequency\n(# infants)")+theme_bw()+xlab("Mean arrival time for each OTU across infants")+ylab("Mean persistence score for each OTU across infants")+theme(axis.text=element_text(size=15))+theme(axis.title=element_text(size=16))+theme(legend.text=element_text(size=12))+theme(legend.title=element_text(size=13))

#figure 2a - taxonomic composition of each priority effects group
colourCount = length(unique(summary_all_OTUs$Class))
ggplot(summary_all_OTUs[summary_all_OTUs$priority_effects!="None" & !is.na(summary_all_OTUs$priority_effects),],aes(priority_effects,fill=Class))+geom_bar()+theme_bw()+xlab("")+scale_fill_manual(values = getPalette(colourCount))+theme(axis.text=element_text(size=15))+theme(axis.title=element_text(size=16))+theme(legend.text=element_text(size=12))+theme(legend.title=element_text(size=13))+ylab("Number of OTUs")

#figure 2b - priority effects distribution of each class
by_class<-data.frame(table(summary_all_OTUs$Class,summary_all_OTUs$priority_effects))
by_class$sum<-rep(table(summary_all_OTUs$Class))
by_class$proportion<-by_class$Freq/by_class$sum
by_class<-by_class[!by_class$Var1%in%by_class[by_class$proportion==1,"Var1"],]
by_class<-by_class[by_class$Var1!="Chloroplast" & by_class$Var1!="unclassified",]
ggplot(by_class[by_class$Var1!="Chloroplast",],aes(Var2,proportion,fill=Var1))+geom_bar(stat="identity")+facet_wrap(~Var1)+theme_bw()+guides(fill=F)+xlab("")+ylab("Proportion of OTUs in each category")+theme(axis.text.x = element_text(angle=45,hjust=1))+theme(axis.text=element_text(size=14))+theme(axis.title=element_text(size=16))+theme(strip.text = element_text(size=16))

#figure 3a - priority effects vs population frequency
ggplot(summary_all_OTUs,aes(num_infants,arrival_persistence_cor,color=arrival_persistence_cor<0))+geom_point()+stat_smooth(method="lm")+guides(color=F)+xlab("Population frequency (# infants)")+ylab("Pearson's correlation between arrival and persistence")+theme_bw()+scale_color_manual(values=c("royalblue","indianred"))+theme(axis.text=element_text(size=15))+theme(axis.title=element_text(size=16))
c

#figure 3b - priority effects vs population frequency
ggplot(summary_all_OTUs[!is.na(summary_all_OTUs$priority_effects),],aes(priority_effects,num_infants,fill=priority_effects))+geom_boxplot()+guides(fill=F)+ylab("Population frequency (# infants)")+theme_bw()+xlab("")+scale_fill_manual(values=c("indianred","grey50","royalblue"))+theme(axis.text=element_text(size=15))+theme(axis.title=element_text(size=16))

#figure 4a - priority effects vs persistence
ggplot(summary_all_OTUs,aes(mean_persistence,arrival_persistence_cor,color=arrival_persistence_cor<0))+geom_point()+stat_smooth(method="lm")+guides(color=F)+xlab("Mean persistence score for each OTU across infants")+ylab("Pearson's correlation between arrival and persistence")+theme_bw()+scale_color_manual(values=c("royalblue","indianred"))+theme(axis.text=element_text(size=15))+theme(axis.title=element_text(size=16))

#figure 4b - priority effects vs persistence
ggplot(summary_all_OTUs[!is.na(summary_all_OTUs$priority_effects),],aes(priority_effects,mean_persistence,fill=priority_effects))+geom_boxplot()+guides(fill=F)+ylab("Mean persistence score for each OTU across infants")+theme_bw()+xlab("")+scale_fill_manual(values=c("indianred","grey50","royalblue"))+theme(axis.text=element_text(size=15))+theme(axis.title=element_text(size=16))

#figure 5a - one-way preemption effects & persistence
ggplot(one_way_preemption,aes(focal_persistence,non_focal_persistence))+geom_point()+theme_bw()+geom_abline(slope=1,linetype="dashed")+xlim(0,1)+ylim(0,1)+xlab("Mean persistence of 'inferior competitor'")+ylab("Mean persistence of 'superior competitor'")+theme(axis.text=element_text(size=15))+theme(axis.title=element_text(size=16))

#figure 5b - one-way preemption effects & population frequency
ggplot(one_way_preemption,aes(focal_occurrence,non_focal_occurrence))+geom_point()+theme_bw()+geom_abline(slope=1,linetype="dashed")+xlab("Population frequency (# infants) of 'inferior competitor'")+ylab("Population frequency (# infants) of 'superior competitor'")+theme(axis.text=element_text(size=15))+theme(axis.title=element_text(size=16))+xlim(0,60)+ylim(0,60)

#figure 6a
grid.arrange(ggplot(summary_all_pairs[summary_all_pairs$preemption=="Arrival time",],aes("",fill=same_phylum))+geom_bar()+theme_bw()+scale_fill_brewer(palette="Set2")+xlab("")+coord_polar("y")+guides(fill=F)+ylab("")+theme(axis.ticks=element_blank())+theme(axis.text=element_text(size=14)),ggplot(summary_all_pairs[summary_all_pairs$preemption=="A or B",],aes("",fill=same_phylum))+geom_bar()+theme_bw()+scale_fill_brewer(palette="Set2")+xlab("")+coord_polar("y")+guides(fill=F)+ylab("")+theme(axis.ticks=element_blank())+theme(axis.text=element_text(size=14)),ggplot(summary_all_pairs[summary_all_pairs$preemption=="Both A and B",],aes("",fill=same_phylum))+geom_bar()+theme_bw()+scale_fill_brewer(palette="Set2")+xlab("")+coord_polar("y")+guides(fill=F)+ylab("")+theme(axis.text=element_text(size=14))+theme(axis.ticks=element_blank()),nrow=1)
grid.arrange(ggplot(summary_all_pairs[summary_all_pairs$preemption=="Arrival time",],aes("",fill=same_class))+geom_bar()+theme_bw()+scale_fill_brewer(palette="Set2")+xlab("")+coord_polar("y")+guides(fill=F)+ylab("")+theme(axis.ticks=element_blank())+theme(axis.text=element_text(size=14)),ggplot(summary_all_pairs[summary_all_pairs$preemption=="A or B",],aes("",fill=same_class))+geom_bar()+theme_bw()+scale_fill_brewer(palette="Set2")+xlab("")+coord_polar("y")+guides(fill=F)+ylab("")+theme(axis.ticks=element_blank())+theme(axis.text=element_text(size=14)),ggplot(summary_all_pairs[summary_all_pairs$preemption=="Both A and B",],aes("",fill=same_class))+geom_bar()+theme_bw()+scale_fill_brewer(palette="Set2")+xlab("")+coord_polar("y")+guides(fill=F)+ylab("")+theme(axis.text=element_text(size=14))+theme(axis.ticks=element_blank()),nrow=1)

#figure 6b
ggplot(summary_all_pairs,aes(avg_co_occurrence,fill=preemption))+geom_density(alpha=0.5)+facet_wrap(~preemption,nrow=3)+theme_bw()+guides(fill=F)+xlab("Average co-occurrence within infants")+theme(axis.text=element_text(size=14))+theme(axis.title=element_text(size=16))+theme(strip.text=element_text(size=15))


