##Arrival & persistence analysis for priority effects review
##Reena Debray

##import libraries
library(ggplot2)
library(readr)
library(RColorBrewer)
library(gridExtra)
library(ape)
library(vegan)
library(phyloseq)
library(DESeq2)
library(dplyr)

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


-----------------------------------------------------------------------------------------------------------------------------------------


###Human gut microbiome data analysis
###data from Guitar et al 2019, Nature Communications

##import OTU tables and OTU taxonomy
#human data: OTU tables ("humangut_otus"), taxonomic annotations ("humangut_taxonomy"), trait annotations ("predicted_trait_data")
humangut_otus <- read_csv("~/Desktop/Microbiomes as food webs/Human gut data/otus.csv")
humangut_taxonomy <- read_csv("~/Desktop/Microbiomes as food webs/Human gut data/SILVA_taxonomy.csv")
predicted_trait_data <- read_csv("~/Desktop/Microbiomes as food webs/Human gut data/predicted_trait_data.csv")


##Make a list of OTUs that are detected in at least 20% of hosts
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


##identify OTUs that are sensitive to microbiome composition
summary_humangut_OTUs<-data.frame(matrix(nrow=0,ncol=13)) #this dataframe contains summary statistics of mean arrival time, mean persistence, and sensitivity to microbiome composition of each OTU
humangut_list<-list() #each entry of this list contains the arrival time & persistence of the OTU in each host
i=1

for (OTU in humangut_OTU_list){
  
  #generate a data frame with the arrival time & persistence of this OTU in each host
  data_per_OTU<-data.frame(matrix(nrow=0,ncol=4))
  for (subjectID in levels(factor(otus$subject))){
    otus_subset<-df_subset(humangut_otus,OTU,subjectID)
    arrival_time<-calculate_arrival_time(otus_subset,OTU)
    persistence<-calculate_persistence(otus_subset,OTU,arrival_time,6) #presence/absence persistence within 6-month window after arrival
    data_per_OTU<-rbind(data_per_OTU,data.frame(OTU,subjectID,arrival_time,persistence))
  }
  
  #calculate summary statistics; find taxonomy
  mean_arrival<-mean(data_per_OTU$arrival_time,na.rm=T)
  mean_persistence<-mean(data_per_OTU$persistence,na.rm=T)
  sd_persistence<-sd(data_per_OTU$persistence,na.rm=T)
  num_hosts<-nrow(data_per_OTU[!is.na(data_per_OTU$persistence),])
  taxonomy<-humangut_taxonomy[humangut_taxonomy$otu==OTU,c(2,4:9)]
  data_per_OTU<-data_per_OTU[!is.na(data_per_OTU$arrival_time) & !is.na(data_per_OTU$persistence),] #arrival time is NA if it never occurred in that host; persistence is NA if no samples were taken within 6 months after arrival
  
  #if there are enough observations, test for dependence on microbiome composition at arrival (t=0) and before arrival (t=-1)
  if (num_hosts>=12 & sd(data_per_OTU$persistence)>0){
  
      ##dependence on microbiome composition at arrival (t=0)
      current_comm_analysis<-data.frame(matrix(nrow=0,ncol=2419))
      colnames(current_comm_analysis)=c("subject","arrival_time","persistence",colnames(otus[5:2420]))
      for (row in seq(1,nrow(data_per_OTU))){
          subject<-as.character(data_per_OTU[row,"subjectID"])
          arrival_time<-data_per_OTU[row,"arrival_time"]
          persistence<-data_per_OTU[row,"persistence"]
  
          #extract community at arrival time in each host
          current_community<-data.frame(otus[otus$subject==subject & otus$t==arrival_time,5:2420])
          current_comm_analysis[row,1:3]<-c(subject,arrival_time,persistence)
          current_comm_analysis[row,4:2419]<-current_community
          }
    
      #remove focal OTU from community
      current_comm_analysis$persistence<-as.numeric(current_comm_analysis$persistence)
      focal<-match(OTU,colnames(current_comm_analysis))
      current_comm_analysis<-current_comm_analysis[,-c(focal)]
  
      #record results of adonis test
      dist<-vegdist(current_comm_analysis[,4:2418], method="bray")
      a<-adonis(dist~persistence,permutations=1000,data=current_comm_analysis)
      persistence_F<-a$aov.tab[1,4]
      persistence_pvalue<-a$aov.tab[1,6]
    
      ##dependence on microbiome composition before arrival (t=-1)
      prior_comm_analysis<-data.frame(matrix(nrow=0,ncol=2419))
      colnames(prior_comm_analysis)=c("subject","arrival_time","persistence",colnames(otus[5:2420]))
      for (row in seq(1,nrow(data_per_OTU))){
          subject<-as.character(data_per_OTU[row,"subjectID"])
          arrival_time<-data_per_OTU[row,"arrival_time"]
          persistence<-data_per_OTU[row,"persistence"]
        
          #in each host, find the last timepoint taken before arrival & extract community
          comm_subset<-otus[otus$subject==subject,]
          comm_subset<-comm_subset[order(comm_subset$t),]
          if (match(arrival_time,comm_subset$t)>1){
              prior_timepoint<-as.numeric(comm_subset[match(arrival_time,comm_subset$t)-1,"t"])
              prior_community<-data.frame(comm_subset[comm_subset$t==prior_timepoint,5:2420])
              prior_comm_analysis[row,1:3]<-c(subject,arrival_time,persistence)
              prior_comm_analysis[row,4:2419]<-prior_community
          }
       }
    
      #remove focal OTU from community
      prior_comm_analysis$persistence<-as.numeric(prior_comm_analysis$persistence)
      focal<-match(OTU,colnames(prior_comm_analysis))
      prior_comm_analysis<-prior_comm_analysis[,-c(focal)]
      prior_comm_analysis<-prior_comm_analysis[!is.na(prior_comm_analysis$arrival_time),]
    
      #record results of adonis test
      dist<-vegdist(prior_comm_analysis[,4:2418], method="bray")
      a<-adonis(dist~persistence,permutations=1000,data=prior_comm_analysis)
      persistence_F_prior<-a$aov.tab[1,4]
      persistence_pvalue_prior<-a$aov.tab[1,6]
  
      #iteratively update data frames
      summary_humangut_OTUs<-rbind(summary_humangut_OTUs,data.frame(OTU,mean_arrival,mean_persistence,sd_persistence,num_hosts,taxonomy,persistence_F,persistence_pvalue,persistence_F_prior,persistence_pvalue_prior))
      humangut_list[[i]]<-data_per_OTU
      i=i+1
     }  
  }



##DESeq analysis of "partner" OTUs whose abundance at (t=0) or before (t=-1) arrival predicts persistence of "focal" OTU

#to be detectably sensitive to microbiome composition, OTUs need to vary a lot in persistence
#here I am considering the set of OTUs with variation in persistence (SD) above the median level-- but this could be adjusted

sd_persist_med<-summary(summary_humangut_OTUs$sd_persistence)["Median"]

tmp<-summary_humangut_OTUs[summary_humangut_OTUs$sd_persistence>=sd_persist_med,]
tmp$persistence_padj<-p.adjust(tmp$persistence_pvalue,method="BH")
humangut_sensitive_OTUs_current<-as.character(tmp[tmp$persistence_padj<0.1,"OTU"])

tmp$persistence_padj_prior<-p.adjust(tmp$persistence_pvalue_prior,method="BH")
humangut_sensitive_OTUs_prior<-as.character(tmp[tmp$persistence_padj_prior<0.1,"OTU"])


##identify partner OTUs within the community at arrival (t=0)
humangut_current_deseq<-data.frame(matrix(nrow=0,ncol=22))

#set up dataframe of community at arrival again (for the small number of sensitive OTUs)
for (focal_OTU in humangut_sensitive_OTUs_current){
  data_per_OTU<-humangut_list[[match(focal_OTU,summary_humangut_OTUs$OTU)]]
  current_comm_analysis<-data.frame(matrix(nrow=0,ncol=2419))
  colnames(current_comm_analysis)=c("subject","arrival_time","persistence",colnames(otus[5:2420]))
  data_per_OTU<-data_per_OTU[!is.na(data_per_OTU$arrival_time) & !is.na(data_per_OTU$persistence),]
  
  for (row in seq(1,nrow(data_per_OTU))){
    subject<-as.character(data_per_OTU[row,"subjectID"])
    arrival_time<-data_per_OTU[row,"arrival_time"]
    persistence<-data_per_OTU[row,"persistence"]
  
    current_community<-data.frame(otus[otus$subject==subject & otus$t==arrival_time,5:2420])
    current_comm_analysis[row,1:3]<-c(subject,arrival_time,persistence)
    current_comm_analysis[row,4:2419]<-current_community
    }

  focal<-match(OTU,colnames(current_comm_analysis))
  current_comm_analysis<-current_comm_analysis[,-c(focal)]
  current_comm_analysis$persistence<-as.numeric(current_comm_analysis$persistence)

  #convert to phyloseq object
  #1: OTU table (convert to numeric, add 1 to everything so DEseq can run)
  phyloseq_otu_table<-t(current_comm_analysis[,4:2418])
  colnames(phyloseq_otu_table)<-current_comm_analysis$subject
  phyloseq_otu_table<-mutate_all(data.frame(phyloseq_otu_table), function(x) as.numeric(as.character(x)))
  phyloseq_otu_table<-mutate_all(data.frame(phyloseq_otu_table), function(x){x+1})
  phyloseq_otu_table<-as.matrix(phyloseq_otu_table)
  phyloseq_otu_table = otu_table(phyloseq_otu_table, taxa_are_rows = TRUE)
  
  #2: taxonomy (remove OTUs that are not present in the OTU table)
  phyloseq_taxonomy<-data.frame(humangut_taxonomy[,-c(1)])
  phyloseq_taxonomy<-phyloseq_taxonomy[!phyloseq_taxonomy$otu%in%setdiff(phyloseq_taxonomy$otu,colnames(current_comm_analysis[,4:2418])),]
  phyloseq_taxonomy<-phyloseq_taxonomy[order(phyloseq_taxonomy$otu),]
  rownames(phyloseq_taxonomy)<-taxa_names(phyloseq_otu_table)
  phyloseq_taxonomy<-tax_table(as.matrix(phyloseq_taxonomy))

  #3: sample info (add a discrete time variable, add rownames)
  phyloseq_samples<-data.frame(current_comm_analysis[,c(1:3)])
  rownames(phyloseq_samples)<-phyloseq_samples[,1]
  phyloseq_samples<-phyloseq_samples[,-c(1)]
  phyloseq_samples<-sample_data(phyloseq_samples)

  #create a phyloseq object
  phyloseq_obj<-phyloseq(phyloseq_otu_table,phyloseq_taxonomy,phyloseq_samples)
  deseq_test = phyloseq_to_deseq2(phyloseq_obj, ~persistence)
  deseq_output<-DESeq(deseq_test, test="Wald", fitType="parametric")
  #extract significantly enriched/depleted taxa
  res = results(deseq_output, cooksCutoff = FALSE)
  sigtab = res[which(res$padj < 0.05), ]
  sigtab = cbind(as(sigtab, "data.frame"), as(tax_table(phyloseq_obj)[rownames(sigtab), ], "matrix"))
  sigtab<-sigtab[order(-sigtab$log2FoldChange),]
  sigtab$x<-seq(1,nrow(sigtab))
  
  #add focal OTU info and bind to master dataframe
  focal_taxonomy<-data.frame(humangut_taxonomy[humangut_taxonomy$otu==focal_OTU,])
  sigtab[,"focal_OTU"]<-focal_taxonomy[,"otu"]
  sigtab[,"focal_Phylum"]<-focal_taxonomy[,"Phylum"]
  sigtab[,"focal_Class"]<-focal_taxonomy[,"Class"]
  sigtab[,"focal_Order"]<-focal_taxonomy[,"Order"]
  sigtab[,"focal_Family"]<-focal_taxonomy[,"Family"]
  sigtab[,"focal_Genus"]<-focal_taxonomy[,"Genus"]
  sigtab[,"focal_Species"]<-focal_taxonomy[,"Species"]
  humangut_current_deseq<-rbind(humangut_current_deseq,sigtab)
}


##identify partner OTUs within the community before arrival (t=-1)                               
humangut_prior_deseq<-data.frame(matrix(nrow=0,ncol=22))
                                 
#set up dataframe of community at arrival again (for the small number of sensitive OTUs)
for (focal_OTU in humangut_sensitive_OTUs_prior){
    data_per_OTU<-humangut_list[[match(focal_OTU,summary_humangut_OTUs$OTU)]]
    prior_comm_analysis<-data.frame(matrix(nrow=12,ncol=2419))
    colnames(prior_comm_analysis)=c("subject","arrival_time","persistence",colnames(otus[5:2420]))
    data_per_OTU<-data_per_OTU[!is.na(data_per_OTU$arrival_time) & !is.na(data_per_OTU$persistence),]
    
    for (row in seq(1,nrow(data_per_OTU))){
      subject<-as.character(data_per_OTU[row,"subjectID"])
      arrival_time<-data_per_OTU[row,"arrival_time"]
      persistence<-data_per_OTU[row,"persistence"]
      comm_subset<-otus[otus$subject==subject,]
      comm_subset<-comm_subset[order(comm_subset$t),]
      if (match(arrival_time,comm_subset$t)>1){
        prior_timepoint<-as.numeric(comm_subset[match(arrival_time,comm_subset$t)-1,"t"])
        prior_community<-data.frame(comm_subset[comm_subset$t==prior_timepoint,5:2420])
        prior_comm_analysis[row,1:3]<-c(subject,arrival_time,persistence)
        prior_comm_analysis[row,4:2419]<-prior_community
       }
    }
  
    prior_comm_analysis$persistence<-as.numeric(prior_comm_analysis$persistence)
    focal<-match(OTU,colnames(prior_comm_analysis))
    prior_comm_analysis<-prior_comm_analysis[,-c(focal)]
    prior_comm_analysis<-prior_comm_analysis[!is.na(prior_comm_analysis$arrival_time),]
  
  #convert to phyloseq object
  #1: OTU table (convert to numeric, add 1 to everything so DEseq can run)
  phyloseq_otu_table<-t(prior_comm_analysis[,4:2418])
  colnames(phyloseq_otu_table)<-prior_comm_analysis$subject
  phyloseq_otu_table<-mutate_all(data.frame(phyloseq_otu_table), function(x) as.numeric(as.character(x)))
  phyloseq_otu_table<-mutate_all(data.frame(phyloseq_otu_table), function(x){x+1})
  phyloseq_otu_table<-as.matrix(phyloseq_otu_table)
  phyloseq_otu_table = otu_table(phyloseq_otu_table, taxa_are_rows = TRUE)
  
  #2: taxonomy (remove OTUs that are not present in the OTU table)
  phyloseq_taxonomy<-data.frame(humangut_taxonomy[,-c(1)])
  phyloseq_taxonomy<-phyloseq_taxonomy[!phyloseq_taxonomy$otu%in%setdiff(phyloseq_taxonomy$otu,colnames(prior_comm_analysis[,4:2418])),]
  phyloseq_taxonomy<-phyloseq_taxonomy[order(phyloseq_taxonomy$otu),]
  rownames(phyloseq_taxonomy)<-taxa_names(phyloseq_otu_table)
  phyloseq_taxonomy<-tax_table(as.matrix(phyloseq_taxonomy))

  #3: sample info (add a discrete time variable, add rownames)
  phyloseq_samples<-data.frame(prior_comm_analysis[,c(1:3)])
  rownames(phyloseq_samples)<-phyloseq_samples[,1]
  phyloseq_samples<-phyloseq_samples[,-c(1)]
  phyloseq_samples<-sample_data(phyloseq_samples)

  #create a phyloseq object
  phyloseq_obj<-phyloseq(phyloseq_otu_table,phyloseq_taxonomy,phyloseq_samples)
  deseq_test = phyloseq_to_deseq2(phyloseq_obj, ~persistence)
  deseq_output<-DESeq(deseq_test, test="Wald", fitType="parametric")
  #extract significantly enriched/depleted taxa
  res = results(deseq_output, cooksCutoff = FALSE)
  sigtab = res[which(res$padj < 0.05), ]
  sigtab = cbind(as(sigtab, "data.frame"), as(tax_table(phyloseq_obj)[rownames(sigtab), ], "matrix"))
  sigtab<-sigtab[order(-sigtab$log2FoldChange),]
  sigtab$x<-seq(1,nrow(sigtab))
  
  #add focal OTU info and bind to master data frame
  focal_taxonomy<-data.frame(humangut_taxonomy[humangut_taxonomy$otu==focal_OTU,])
  sigtab[,"focal_OTU"]<-focal_taxonomy[,"otu"]
  sigtab[,"focal_Phylum"]<-focal_taxonomy[,"Phylum"]
  sigtab[,"focal_Class"]<-focal_taxonomy[,"Class"]
  sigtab[,"focal_Order"]<-focal_taxonomy[,"Order"]
  sigtab[,"focal_Family"]<-focal_taxonomy[,"Family"]
  sigtab[,"focal_Genus"]<-focal_taxonomy[,"Genus"]
  sigtab[,"focal_Species"]<-focal_taxonomy[,"Species"]
  humangut_prior_deseq<-rbind(humangut_prior_deseq,sigtab)
}


##Compare relatedness of DESeq-identified pairs to 1000 iterations of random pairs in the human gut microbiome
#How many DESeq-identified pairs come from the same family? (t=-1)
humangut_perm_current<-c()
all_pairs<-combn(as.character(summary_humangut_OTUs$OTU),2) #all possible pairs of OTUs that passed our initial filtering criteria
for (i in seq(1,1000)){
  print(i)
  null_relatedness_data<-data.frame(matrix(nrow=nrow(humangut_current_deseq),ncol=5))
  colnames(null_relatedness_data)<-c("OTU_A","Family_A","OTU_B","Family_B","relatedness")
  sampled_pairs<-all_pairs[,sample(seq(1,ncol(all_pairs)),nrow(humangut_current_deseq))] #randomly sample the same number of pairs as in the observe data

  for (j in seq(1,nrow(humangut_current_deseq))){
    OTU_A<-sampled_pairs[1,j]
    family_A<-as.character(humangut_taxonomy[humangut_taxonomy$otu==OTU_A,7])
    OTU_B<-sampled_pairs[2,j]
    family_B<-as.character(humangut_taxonomy[humangut_taxonomy$otu==OTU_B,7])

    if (family_A!="unclassified" & family_B!="unclassified"){
    relatedness<-as.numeric(family_A==family_B)
    null_relatedness_data[j,]<-c(OTU_A,family_A,OTU_B,family_B,relatedness)
    }
  }
  null_relatedness_data$relatedness<-as.numeric(null_relatedness_data$relatedness)
  humangut_perm_current<-c(humangut_perm_current,mean(null_relatedness_data$relatedness,na.rm=T))
}

obs_means_neg<-nrow(humangut_current_deseq[humangut_current_deseq$log2FoldChange<0 & humangut_current_deseq$focal_Family!="unclassified" & humangut_current_deseq$Family!="unclassified" & humangut_current_deseq$Family==humangut_current_deseq$focal_Family,])/nrow(humangut_current_deseq[humangut_current_deseq$log2FoldChange<0 & humangut_current_deseq$focal_Family!="unclassified" & humangut_current_deseq$Family!="unclassified",])
print(length(humangut_perm_current[humangut_perm_current<=obs_means_neg])/1000)
print(length(humangut_perm_current[humangut_perm_current>=obs_means_neg])/1000)

                                 
obs_means_pos<-nrow(humangut_current_deseq[humangut_current_deseq$log2FoldChange>0 & humangut_current_deseq$focal_Family!="unclassified" & humangut_current_deseq$Family!="unclassified" & humangut_current_deseq$Family==humangut_current_deseq$focal_Family,])/nrow(humangut_current_deseq[humangut_current_deseq$log2FoldChange>0 & humangut_current_deseq$focal_Family!="unclassified" & humangut_current_deseq$Family!="unclassified",])
print(length(humangut_perm_current[humangut_perm_current<=obs_means_pos])/1000)
print(length(humangut_perm_current[humangut_perm_current>=obs_means_pos])/1000)

#positive pairs from arrival (t=0) are slightly MORE closely related than expected by chance (p=0.032)                               
                                 
#How many DESeq-identified pairs come from the same family? (t=-1)
humangut_perm_prior<-c()
all_pairs<-combn(as.character(summary_humangut_OTUs$OTU),2) #all possible pairs of OTUs that passed our initial filtering criteria
for (i in seq(1,1000)){
  print(i)
  null_relatedness_data<-data.frame(matrix(nrow=nrow(humangut_prior_deseq),ncol=5))
  colnames(null_relatedness_data)<-c("OTU_A","Family_A","OTU_B","Family_B","relatedness")
  sampled_pairs<-all_pairs[,sample(seq(1,ncol(all_pairs)),nrow(humangut_prior_deseq))] #randomly sample the same number of pairs as in the observe data

  for (j in seq(1,nrow(humangut_prior_deseq))){
    OTU_A<-sampled_pairs[1,j]
    family_A<-as.character(humangut_taxonomy[humangut_taxonomy$otu==OTU_A,7])
    OTU_B<-sampled_pairs[2,j]
    family_B<-as.character(humangut_taxonomy[humangut_taxonomy$otu==OTU_B,7])

    if (family_A!="unclassified" & family_B!="unclassified"){
    relatedness<-as.numeric(family_A==family_B)
    null_relatedness_data[j,]<-c(OTU_A,family_A,OTU_B,family_B,relatedness)
    }
  }
  null_relatedness_data$relatedness<-as.numeric(null_relatedness_data$relatedness)
  humangut_perm_prior<-c(humangut_perm_prior,mean(null_relatedness_data$relatedness,na.rm=T))
}
                                 
obs_means_neg<-nrow(humangut_prior_deseq[humangut_prior_deseq$log2FoldChange<0 & humangut_prior_deseq$focal_Family!="unclassified" & humangut_prior_deseq$Family!="unclassified" & humangut_prior_deseq$Family==humangut_prior_deseq$focal_Family,])/nrow(humangut_prior_deseq[humangut_prior_deseq$log2FoldChange<0 & humangut_prior_deseq$focal_Family!="unclassified" & humangut_prior_deseq$Family!="unclassified",])
print(length(humangut_perm_prior[humangut_perm_prior<=obs_means_neg])/1000)
print(length(humangut_perm_prior[humangut_perm_prior>=obs_means_neg])/1000)
                                 
#negative pairs from the time before arrival (t=-1) are MORE closely related than expected by chance (p=0.001)
                                 
obs_means_pos<-nrow(humangut_prior_deseq[humangut_prior_deseq$log2FoldChange>0 & humangut_prior_deseq$focal_Family!="unclassified" & humangut_prior_deseq$Family!="unclassified" & humangut_prior_deseq$Family==humangut_prior_deseq$focal_Family,])/nrow(humangut_prior_deseq[humangut_prior_deseq$log2FoldChange>0 & humangut_prior_deseq$focal_Family!="unclassified" & humangut_prior_deseq$Family!="unclassified",])
print(length(humangut_perm_prior[humangut_perm_prior<=obs_means_pos])/1000)
print(length(humangut_perm_prior[humangut_perm_prior>=obs_means_pos])/1000)
                             
                                                                
                                 
-----------------------------------------------------------------------------------------------------------------------------------------                               
                  
                                 
###Mouse gut microbiome data analysis
###Data from Kozich et al 2013, Applied and Environmental Microbiology

##Import OTU tables and OTU taxonomy
murine_otus<-read_csv("~/Desktop/Microbiomes as food webs/Mouse gut data/murine_OTUs.csv")
murine_taxonomy <- read_excel("~/Desktop/Microbiomes as food webs/Mouse gut data/murine_taxonomy.xlsx")
                                                           
##Make a list of OTUs that are detected in at least 20% of hosts
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
                                 
##identify OTUs that are sensitive to microbiome composition
summary_murine_OTUs<-data.frame(matrix(nrow=0,ncol=13)) #this dataframe contains summary statistics of mean arrival time, mean persistence, and sensitivity to microbiome composition of each OTU
murine_list<-list() #each entry of this list contains the arrival time & persistence of the OTU in each host
i=1

for (OTU in murine_OTU_list){
  
  #generate a data frame with the arrival time & persistence of this OTU in each host
  data_per_OTU<-data.frame(matrix(nrow=0,ncol=4))
  for (subjectID in levels(factor(murine_otus$subject))){
    otus_subset<-df_subset(murine_otus,OTU,subjectID)
    arrival_time<-calculate_arrival_time(otus_subset,OTU)
    persistence<-calculate_persistence(otus_subset,OTU,arrival_time,30) #1 month rather than 6 months?
    data_per_OTU<-rbind(data_per_OTU,data.frame(OTU,subjectID,arrival_time,persistence))
  }
  
  #calculate summary statistics; find taxonomy
  mean_arrival<-mean(data_per_OTU$arrival_time,na.rm=T)
  mean_persistence<-mean(data_per_OTU$persistence,na.rm=T)
  sd_persistence<-sd(data_per_OTU$persistence,na.rm=T)
  num_hosts<-nrow(data_per_OTU[!is.na(data_per_OTU$persistence),])
  taxonomy<-murine_taxonomy[murine_taxonomy$OTU==OTU,c(5,7,9,11,13)]
  data_per_OTU<-data_per_OTU[!is.na(data_per_OTU$arrival_time) & !is.na(data_per_OTU$persistence),]

  #if there are enough observations, test for dependence on microbiome composition at arrival (t=0) and before arrival (t=-1)
  if (num_hosts>=3 & sd(data_per_OTU$persistence)>0){
    
      ##dependence on microbiome composition at arrival
      current_comm_analysis<-data.frame(matrix(nrow=0,ncol=3098))
      colnames(current_comm_analysis)=c("subject","arrival_time","persistence",colnames(murine_otus[5:3099]))
      for (row in seq(1,nrow(data_per_OTU))){
         subject<-as.character(data_per_OTU[row,"subjectID"])
         arrival_time<-data_per_OTU[row,"arrival_time"]
         persistence<-data_per_OTU[row,"persistence"]
  
          #extract community at arrival time in each host
          current_community<-data.frame(murine_otus[murine_otus$subject==subject & murine_otus$t==arrival_time,5:3099])
          current_comm_analysis[row,1:3]<-c(subject,arrival_time,persistence)
          current_comm_analysis[row,4:3098]<-current_community
          }
    
       #remove focal OTU from community
       current_comm_analysis$persistence<-as.numeric(current_comm_analysis$persistence)
       focal<-match(OTU,colnames(current_comm_analysis))
       current_comm_analysis<-current_comm_analysis[,-c(focal)]
  
       #record results of adonis test
       dist<-vegdist(current_comm_analysis[,3:3097], method="bray")
       a<-adonis(dist~persistence,permutations=1000,data=current_comm_analysis)
       persistence_F<-a$aov.tab[1,4]
       persistence_pvalue<-a$aov.tab[1,6]
    
       ##dependence on microbiome composition before arrival
       prior_comm_analysis<-data.frame(matrix(nrow=0,ncol=3098))
       colnames(prior_comm_analysis)=c("subject","arrival_time","persistence",colnames(murine_otus[5:3099]))
       for (row in seq(1,nrow(data_per_OTU))){
          subject<-as.character(data_per_OTU[row,"subjectID"])
          arrival_time<-data_per_OTU[row,"arrival_time"]
          persistence<-data_per_OTU[row,"persistence"]
        
          #in each host, find the last timepoint taken before arrival & extract community
          comm_subset<-murine_otus[murine_otus$subject==subject,]
          comm_subset<-comm_subset[order(comm_subset$t),]
          if (match(arrival_time,comm_subset$t)>1){
              prior_timepoint<-as.numeric(comm_subset[match(arrival_time,comm_subset$t)-1,"t"])
              prior_community<-data.frame(comm_subset[comm_subset$t==prior_timepoint,5:3099])
              prior_comm_analysis[row,1:3]<-c(subject,arrival_time,persistence)
              prior_comm_analysis[row,4:3098]<-prior_community
          }
      }
    
      #remove focal OTU from community
      prior_comm_analysis$persistence<-as.numeric(prior_comm_analysis$persistence)
      focal<-match(OTU,colnames(prior_comm_analysis))
      prior_comm_analysis<-prior_comm_analysis[,-c(focal)]
      prior_comm_analysis<-prior_comm_analysis[!is.na(prior_comm_analysis$arrival_time),]
  
      #record results of adonis test
      if (nrow(prior_comm_analysis[prior_comm_analysis$persistence!=1 & prior_comm_analysis$persistence!=0,])>2){
          dist<-vegdist(prior_comm_analysis[,3:3097], method="bray")
          a<-adonis(dist~persistence,permutations=1000,data=prior_comm_analysis)
          persistence_F_prior<-a$aov.tab[1,4]
          persistence_pvalue_prior<-a$aov.tab[1,6]
        }
      else {persistence_F_prior<-NA; persistence_pvalue_prior<-NA}
    
      summary_murine_OTUs<-rbind(summary_murine_OTUs,data.frame(OTU,mean_arrival,mean_persistence,sd_persistence,num_hosts,taxonomy,persistence_F,persistence_pvalue,persistence_F_prior,persistence_pvalue_prior))
      murine_list[[i]]<-data_per_OTU
      i=i+1
  }
 
}


##DESeq analysis of "partner" OTUs whose abundance at (t=0) or before (t=-1) arrival predicts persistence of "focal" OTU

#to be detectably sensitive to microbiome composition, OTUs need to vary a lot in persistence
#here I am considering the set of OTUs with variation in persistence (SD) above the median level-- but this could be adjusted

sd_persist_med<-summary(summary_murine_OTUs$sd_persistence)["Median"]

tmp<-summary_murine_OTUs[summary_murine_OTUs$sd_persistence>=sd_persist_med,]
tmp$persistence_padj<-p.adjust(tmp$persistence_pvalue,method="BH")
murine_sensitive_OTUs_current<-as.character(tmp[tmp$persistence_padj<0.1,"OTU"])

tmp$persistence_padj_prior<-p.adjust(tmp$persistence_pvalue_prior,method="BH")
murine_sensitive_OTUs_prior<-as.character(tmp[tmp$persistence_padj_prior<0.1,"OTU"])

#None of these pass multiple testing correction (likely due to small sample size + low variation in persistence of many OTUs).
#Let's use the set that are nominally significant for now-- to see if patterns are the same as the other two datasets
                                 
murine_sensitive_OTUs_current<-summary_murine_OTUs[summary_murine_OTUs$persistence_pvalue<0.05 & !is.na(summary_murine_OTUs$persistence_pvalue),"OTU"]
murine_sensitive_OTUs_prior<-summary_murine_OTUs[summary_murine_OTUs$persistence_pvalue_prior<0.05 & !is.na(summary_murine_OTUs$persistence_pvalue_prior),"OTU"]
                                 
                                 
##identify partner OTUs within the community at arrival (t=0)
murine_current_deseq<-data.frame(matrix(nrow=0,ncol=22))

#set up dataframe of community at arrival again (for the small number of sensitive OTUs)
for (focal_OTU in murine_sensitive_OTUs_current){
  data_per_OTU<-murine_list[[match(focal_OTU,summary_murine_OTUs$OTU)]]
  current_comm_analysis<-data.frame(matrix(nrow=0,ncol=3098))
  colnames(current_comm_analysis)=c("subject","arrival_time","persistence",colnames(murine_otus[5:3099]))
  data_per_OTU<-data_per_OTU[!is.na(data_per_OTU$arrival_time) & !is.na(data_per_OTU$persistence),]
  
  for (row in seq(1,nrow(data_per_OTU))){
    subject<-as.character(data_per_OTU[row,"subjectID"])
    arrival_time<-data_per_OTU[row,"arrival_time"]
    persistence<-data_per_OTU[row,"persistence"]
    current_community<-data.frame(murine_otus[murine_otus$subject==subject & murine_otus$t==arrival_time,5:3099])
    current_comm_analysis[row,1:3]<-c(subject,arrival_time,persistence)
    current_comm_analysis[row,4:3098]<-current_community
    }

  focal<-match(OTU,colnames(current_comm_analysis))
  current_comm_analysis<-current_comm_analysis[,-c(focal)]
  current_comm_analysis$persistence<-as.numeric(current_comm_analysis$persistence)
  
  #convert to phyloseq object
  #1: OTU table (convert to numeric, add 1 to everything so DEseq can run)
  phyloseq_otu_table<-t(current_comm_analysis[,4:3097])
  colnames(phyloseq_otu_table)<-current_comm_analysis$subject
  phyloseq_otu_table<-mutate_all(data.frame(phyloseq_otu_table), function(x) as.numeric(as.character(x)))
  phyloseq_otu_table<-mutate_all(data.frame(phyloseq_otu_table), function(x){x+1})
  phyloseq_otu_table<-as.matrix(phyloseq_otu_table)
  phyloseq_otu_table = otu_table(phyloseq_otu_table, taxa_are_rows = TRUE)
  
  #2: taxonomy (remove OTUs that are not present in the OTU table)
  phyloseq_taxonomy<-data.frame(murine_taxonomy[,seq(1,13,2)])
  phyloseq_taxonomy<-phyloseq_taxonomy[!phyloseq_taxonomy$OTU%in%setdiff(phyloseq_taxonomy$OTU,colnames(current_comm_analysis[,4:3097])),]
  phyloseq_taxonomy<-phyloseq_taxonomy[order(phyloseq_taxonomy$OTU),]
  rownames(phyloseq_taxonomy)<-taxa_names(phyloseq_otu_table)
  phyloseq_taxonomy<-tax_table(as.matrix(phyloseq_taxonomy))

  #3: sample info (add rownames)
  phyloseq_samples<-data.frame(current_comm_analysis[,c(1:3)])
  rownames(phyloseq_samples)<-phyloseq_samples[,1]
  phyloseq_samples<-phyloseq_samples[,-c(1)]
  phyloseq_samples<-sample_data(phyloseq_samples)
                                                            
  #create a phyloseq object
  phyloseq_obj<-phyloseq(phyloseq_otu_table,phyloseq_taxonomy,phyloseq_samples)
  deseq_test = phyloseq_to_deseq2(phyloseq_obj, ~persistence)
  deseq_output<-DESeq(deseq_test, test="Wald", fitType="mean")
  #extract significantly enriched/depleted taxa
  res = results(deseq_output, cooksCutoff = FALSE)
  sigtab = res[which(res$padj < 0.05), ]
  if (nrow(sigtab)>1){
      sigtab = cbind(as(sigtab, "data.frame"), as(tax_table(phyloseq_obj)[rownames(sigtab), ], "matrix"))
      sigtab<-sigtab[order(-sigtab$log2FoldChange),]
      sigtab$x<-seq(1,nrow(sigtab))
  
      #add focal OTU info and bind to master dataframe
      taxonomy<-as.character(murine_taxonomy[murine_taxonomy$OTU==focal_OTU,])
      sigtab[,"focal_OTU"]<-focal_OTU
      sigtab[,"focal_Phylum"]<-taxonomy[5]
      sigtab[,"focal_Class"]<-taxonomy[7]
      sigtab[,"focal_Order"]<-taxonomy[9]
      sigtab[,"focal_Family"]<-taxonomy[11]
      sigtab[,"focal_Genus"]<-taxonomy[13]
      murine_current_deseq<-rbind(murine_current_deseq,sigtab)
    }
}                                        
                                 
                                 
##identify partner OTUs within the community before arrival (t=-1)
murine_prior_deseq<-data.frame(matrix(nrow=0,ncol=22))

#set up dataframe of community at arrival again (for the small number of sensitive OTUs)
for (focal_OTU in murine_sensitive_OTUs_prior){
  data_per_OTU<-murine_list[[match(focal_OTU,summary_murine_OTUs$OTU)]]
  prior_comm_analysis<-data.frame(matrix(nrow=0,ncol=3098))
  colnames(prior_comm_analysis)=c("subject","arrival_time","persistence",colnames(murine_otus[5:3099]))
  data_per_OTU<-data_per_OTU[!is.na(data_per_OTU$arrival_time) & !is.na(data_per_OTU$persistence),]
  
  for (row in seq(1,nrow(data_per_OTU))){
    subject<-as.character(data_per_OTU[row,"subjectID"])
    arrival_time<-data_per_OTU[row,"arrival_time"]
    persistence<-data_per_OTU[row,"persistence"]
    
    
    #in each host, find the last timepoint taken before arrival & extract community
    comm_subset<-murine_otus[murine_otus$subject==subject,]
    comm_subset<-comm_subset[order(comm_subset$t),]
    if (match(arrival_time,comm_subset$t)>1){
        prior_timepoint<-as.numeric(comm_subset[match(arrival_time,comm_subset$t)-1,"t"])
        prior_community<-data.frame(comm_subset[comm_subset$t==prior_timepoint,5:3099])
        prior_comm_analysis[row,1:3]<-c(subject,arrival_time,persistence)
        prior_comm_analysis[row,4:3098]<-prior_community
        }
     }
    
     #remove focal OTU from community
     prior_comm_analysis$persistence<-as.numeric(prior_comm_analysis$persistence)
     focal<-match(OTU,colnames(prior_comm_analysis))
     prior_comm_analysis<-prior_comm_analysis[,-c(focal)]
     prior_comm_analysis<-prior_comm_analysis[!is.na(prior_comm_analysis$arrival_time),]
  
  #convert to phyloseq object
  #1: OTU table (convert to numeric, add 1 to everything so DEseq can run)
  phyloseq_otu_table<-t(prior_comm_analysis[,4:3097])
  colnames(phyloseq_otu_table)<-prior_comm_analysis$subject
  phyloseq_otu_table<-mutate_all(data.frame(phyloseq_otu_table), function(x) as.numeric(as.character(x)))
  phyloseq_otu_table<-mutate_all(data.frame(phyloseq_otu_table), function(x){x+1})
  phyloseq_otu_table<-as.matrix(phyloseq_otu_table)
  phyloseq_otu_table = otu_table(phyloseq_otu_table, taxa_are_rows = TRUE)
  
  #2: taxonomy (remove OTUs that are not present in the OTU table)
  phyloseq_taxonomy<-data.frame(murine_taxonomy[,seq(1,13,2)])
  phyloseq_taxonomy<-phyloseq_taxonomy[!phyloseq_taxonomy$OTU%in%setdiff(phyloseq_taxonomy$OTU,colnames(prior_comm_analysis[,4:3097])),]
  phyloseq_taxonomy<-phyloseq_taxonomy[order(phyloseq_taxonomy$OTU),]
  rownames(phyloseq_taxonomy)<-taxa_names(phyloseq_otu_table)
  phyloseq_taxonomy<-tax_table(as.matrix(phyloseq_taxonomy))

  #3: sample info (add rownames)
  phyloseq_samples<-data.frame(prior_comm_analysis[,c(1:3)])
  rownames(phyloseq_samples)<-phyloseq_samples[,1]
  phyloseq_samples<-phyloseq_samples[,-c(1)]
  phyloseq_samples<-sample_data(phyloseq_samples)
                                                            
  #create a phyloseq object
  phyloseq_obj<-phyloseq(phyloseq_otu_table,phyloseq_taxonomy,phyloseq_samples)
  deseq_test = phyloseq_to_deseq2(phyloseq_obj, ~persistence)
  deseq_output<-DESeq(deseq_test, test="Wald", fitType="mean")
  #extract significantly enriched/depleted taxa
  res = results(deseq_output, cooksCutoff = FALSE)
  sigtab = res[which(res$padj < 0.05), ]
  if (nrow(sigtab)>=1){
      sigtab = cbind(as(sigtab, "data.frame"), as(tax_table(phyloseq_obj)[rownames(sigtab), ], "matrix"))
      sigtab<-sigtab[order(-sigtab$log2FoldChange),]
      sigtab$x<-seq(1,nrow(sigtab))
  
      #add focal OTU info and bind to master dataframe
      taxonomy<-as.character(murine_taxonomy[murine_taxonomy$OTU==focal_OTU,])
      sigtab[,"focal_OTU"]<-focal_OTU
      sigtab[,"focal_Phylum"]<-taxonomy[5]
      sigtab[,"focal_Class"]<-taxonomy[7]
      sigtab[,"focal_Order"]<-taxonomy[9]
      sigtab[,"focal_Family"]<-taxonomy[11]
      sigtab[,"focal_Genus"]<-taxonomy[13]
      murine_prior_deseq<-rbind(murine_prior_deseq,sigtab)
    }
}                                 
                                 
                                 
                                 
##Compare relatedness of DESeq-identified pairs to 1000 iterations of random pairs in the mouse gut microbiome
#How many DESeq-identified pairs come from the same family? (t=0)
murine_perm_current<-c()
all_pairs<-combn(as.character(summary_murine_OTUs$OTU),2) #all possible pairs of OTUs that passed our initial filtering criteria
for (i in seq(1,1000)){
  print(i)
  null_relatedness_data<-data.frame(matrix(nrow=nrow(murine_current_deseq),ncol=5))
  colnames(null_relatedness_data)<-c("OTU_A","Family_A","OTU_B","Family_B","relatedness")
  sampled_pairs<-all_pairs[,sample(seq(1,ncol(all_pairs)),nrow(murine_current_deseq))] #randomly sample the same number of pairs as in the observed data

  for (j in seq(1,nrow(murine_current_deseq))){
    OTU_A<-sampled_pairs[1,j]
    family_A<-as.character(murine_taxonomy[murine_taxonomy$OTU==OTU_A,"Family"])
    OTU_B<-sampled_pairs[2,j]
    family_B<-as.character(murine_taxonomy[murine_taxonomy$OTU==OTU_B,"Family"])

    if (substr(family_A,nchar(family_A)-11,nchar(family_A))!="unclassified" & substr(family_B,nchar(family_B)-11,nchar(family_B))!="unclassified"){
    relatedness<-as.numeric(family_A==family_B)
    null_relatedness_data[j,]<-c(OTU_A,family_A,OTU_B,family_B,relatedness)
    }
  }
  null_relatedness_data$relatedness<-as.numeric(null_relatedness_data$relatedness)
  murine_perm_current<-c(murine_perm_current,mean(null_relatedness_data$relatedness,na.rm=T))
}

obs_means_neg<-nrow(murine_current_deseq[murine_current_deseq$log2FoldChange<0 & substr(murine_current_deseq$focal_Family,nchar(murine_current_deseq$focal_Family)-11,nchar(murine_current_deseq$focal_Family))!="unclassified" & substr(murine_current_deseq$Family,nchar(murine_current_deseq$Family)-11,nchar(murine_current_deseq$Family))!="unclassified"  & murine_current_deseq$focal_Family==murine_current_deseq$Family,])/nrow(murine_current_deseq[murine_current_deseq$log2FoldChange<0 & substr(murine_current_deseq$focal_Family,nchar(murine_current_deseq$focal_Family)-11,nchar(murine_current_deseq$focal_Family))!="unclassified" & substr(murine_current_deseq$Family,nchar(murine_current_deseq$Family)-11,nchar(murine_current_deseq$Family))!="unclassified",])
print(length(murine_perm_current[murine_perm_current<=obs_means_neg])/1000)
print(length(murine_perm_current[murine_perm_current>=obs_means_neg])/1000)

                                 
obs_means_pos<-nrow(murine_current_deseq[murine_current_deseq$log2FoldChange>0 & substr(murine_current_deseq$focal_Family,nchar(murine_current_deseq$focal_Family)-11,nchar(murine_current_deseq$focal_Family))!="unclassified" & substr(murine_current_deseq$Family,nchar(murine_current_deseq$Family)-11,nchar(murine_current_deseq$Family))!="unclassified"  & murine_current_deseq$focal_Family==murine_current_deseq$Family,])/nrow(murine_current_deseq[murine_current_deseq$log2FoldChange>0 & substr(murine_current_deseq$focal_Family,nchar(murine_current_deseq$focal_Family)-11,nchar(murine_current_deseq$focal_Family))!="unclassified" & substr(murine_current_deseq$Family,nchar(murine_current_deseq$Family)-11,nchar(murine_current_deseq$Family))!="unclassified",])
print(length(murine_perm_current[murine_perm_current<=obs_means_pos])/1000)
print(length(murine_perm_current[murine_perm_current>=obs_means_pos])/1000)

                                 
#How many DESeq-identified pairs come from the same family? (t=-1)
murine_perm_prior<-c()
all_pairs<-combn(as.character(summary_murine_OTUs$OTU),2) #all possible pairs of OTUs that passed our initial filtering criteria
for (i in seq(1,1000)){
  print(i)
  null_relatedness_data<-data.frame(matrix(nrow=nrow(murine_prior_deseq),ncol=5))
  colnames(null_relatedness_data)<-c("OTU_A","Family_A","OTU_B","Family_B","relatedness")
  sampled_pairs<-all_pairs[,sample(seq(1,ncol(all_pairs)),nrow(murine_prior_deseq))] #randomly sample the same number of pairs as in the observed data

  for (j in seq(1,nrow(murine_prior_deseq))){
    OTU_A<-sampled_pairs[1,j]
    family_A<-as.character(murine_taxonomy[murine_taxonomy$OTU==OTU_A,"Family"])
    OTU_B<-sampled_pairs[2,j]
    family_B<-as.character(murine_taxonomy[murine_taxonomy$OTU==OTU_B,"Family"])

    if (substr(family_A,nchar(family_A)-11,nchar(family_A))!="unclassified" & substr(family_B,nchar(family_B)-11,nchar(family_B))!="unclassified"){
    relatedness<-as.numeric(family_A==family_B)
    null_relatedness_data[j,]<-c(OTU_A,family_A,OTU_B,family_B,relatedness)
    }
  }
  null_relatedness_data$relatedness<-as.numeric(null_relatedness_data$relatedness)
  murine_perm_prior<-c(murine_perm_prior,mean(null_relatedness_data$relatedness,na.rm=T))
}

obs_means_neg<-nrow(murine_prior_deseq[murine_prior_deseq$log2FoldChange<0 & substr(murine_prior_deseq$focal_Family,nchar(murine_prior_deseq$focal_Family)-11,nchar(murine_prior_deseq$focal_Family))!="unclassified" & substr(murine_prior_deseq$Family,nchar(murine_prior_deseq$Family)-11,nchar(murine_prior_deseq$Family))!="unclassified"  & murine_prior_deseq$focal_Family==murine_prior_deseq$Family,])/nrow(murine_prior_deseq[murine_prior_deseq$log2FoldChange<0 & substr(murine_prior_deseq$focal_Family,nchar(murine_prior_deseq$focal_Family)-11,nchar(murine_prior_deseq$focal_Family))!="unclassified" & substr(murine_prior_deseq$Family,nchar(murine_prior_deseq$Family)-11,nchar(murine_prior_deseq$Family))!="unclassified",])
print(length(murine_perm_prior[murine_perm_prior<=obs_means_neg])/1000)
print(length(murine_perm_prior[murine_perm_prior>=obs_means_neg])/1000)

                                 
obs_means_pos<-nrow(murine_prior_deseq[murine_prior_deseq$log2FoldChange>0 & substr(murine_prior_deseq$focal_Family,nchar(murine_prior_deseq$focal_Family)-11,nchar(murine_prior_deseq$focal_Family))!="unclassified" & substr(murine_prior_deseq$Family,nchar(murine_prior_deseq$Family)-11,nchar(murine_prior_deseq$Family))!="unclassified"  & murine_prior_deseq$focal_Family==murine_prior_deseq$Family,])/nrow(murine_prior_deseq[murine_prior_deseq$log2FoldChange>0 & substr(murine_prior_deseq$focal_Family,nchar(murine_prior_deseq$focal_Family)-11,nchar(murine_prior_deseq$focal_Family))!="unclassified" & substr(murine_prior_deseq$Family,nchar(murine_prior_deseq$Family)-11,nchar(murine_prior_deseq$Family))!="unclassified",])
print(length(murine_perm_prior[murine_perm_prior<=obs_means_pos])/1000)
print(length(murine_perm_prior[murine_perm_prior>=obs_means_pos])/1000)

-----------------------------------------------------------------------------------------------------------------------------------------                                  
   
                                 
###Cow rumen microbiome
###Data from Furman et al 2020, Nature Communications
                                                               
##Import OTU tables and OTU taxonomy
rumen_otus <- read_csv("~/Desktop/Microbiomes as food webs/Rumen gut data/rumen_OTUs.csv")
rumen_taxonomy <- read_csv("~/Desktop/Microbiomes as food webs/Rumen gut data/Rumen gut taxonomy.csv")

                                 
##Make a list of OTUs that are detected in at least 20% of hosts
rumen_OTU_list<-c()
for (OTU in colnames(rumen_otus[4:2547])){ #loop through all OTUs detected at any point in the study (2544)
  OTU_all_samples<-rumen_otus[,c("subject",OTU)] #subset to this specific OTU
  colnames(OTU_all_samples)<-c("subject","abundance")
  OTU_all_indiv<-OTU_all_samples[OTU_all_samples$abundance>0,] #subset to samples in which this OTU had non-zero abundance
  if (length(unique(OTU_all_indiv$subject))>=9){rumen_OTU_list<-c(rumen_OTU_list,OTU)} #count the number of unique individuals in which this OTU had non-zero abundance for any amount of time
  #consider this OTU only if it is present in 20% or more of individuals (9)
  #Final list includes 2544 OTUs (because the original study had already subsetted to OTUs present in 80% of individuals)
  }
                       

##Identify OTUs that are sensitive to microbiome composition

summary_rumen_OTUs<-data.frame(matrix(nrow=0,ncol=13)) #this dataframe contains summary statistics of mean arrival time, mean persistence, and sensitivity to microbiome composition of each OTU
rumen_list<-list() #each entry of this list contains the arrival time & persistence of the OTU in each host
i=1

for (OTU in rumen_OTU_list){
  data_per_OTU<-data.frame(matrix(nrow=0,ncol=4))
  
  #generate a data frame with the arrival time & persistence of this OTU in each host
  for (subjectID in levels(factor(rumen_otus$subject))){
    otus_subset<-df_subset(rumen_otus,OTU,subjectID)
    arrival_time<-calculate_arrival_time(otus_subset,OTU)
    persistence<-calculate_persistence(otus_subset,OTU,arrival_time,180)  #presence/absence persistence within 6-month window (180 days) after arrival
    data_per_OTU<-rbind(data_per_OTU,data.frame(OTU,subjectID,arrival_time,persistence))
  }
  
  #calculate summary statistics; find taxonomy
  mean_arrival<-mean(data_per_OTU$arrival_time,na.rm=T)
  mean_persistence<-mean(data_per_OTU$persistence,na.rm=T)
  sd_persistence<-sd(data_per_OTU$persistence,na.rm=T)
  num_hosts<-nrow(data_per_OTU[!is.na(data_per_OTU$persistence),])
  taxonomy<-rumen_taxonomy[rumen_taxonomy$X1==OTU,3:8]
  phylum<-substr(taxonomy[1],4,nchar(taxonomy[1]))
  class<-substr(taxonomy[2],4,nchar(taxonomy[2]))
  order<-substr(taxonomy[3],4,nchar(taxonomy[3]))
  family<-substr(taxonomy[4],4,nchar(taxonomy[4]))
  genus<-substr(taxonomy[5],4,nchar(taxonomy[5]))
  species<-substr(taxonomy[6],4,nchar(taxonomy[6]))
  data_per_OTU<-data_per_OTU[!is.na(data_per_OTU$arrival_time) & !is.na(data_per_OTU$persistence),]  #arrival time is NA if it never occurred in that host; persistence is NA if no samples were taken within 6 months after arrival
  
  #Test for dependence on microbiome composition at arrival (t=0) and before arrival (t=-1)
    if (num_hosts>=9 & sd(data_per_OTU$persistence)>0){
      
      ##dependence on microbiome composition at arrival (t=0)
      current_comm_analysis<-data.frame(matrix(nrow=0,ncol=2547))
      colnames(current_comm_analysis)=c("subject","arrival_time","persistence",colnames(rumen_otus[4:2547]))
      for (row in seq(1,nrow(data_per_OTU))){
          subject<-as.character(data_per_OTU[row,"subjectID"])
          arrival_time<-data_per_OTU[row,"arrival_time"]
          persistence<-data_per_OTU[row,"persistence"]
  
          #extract community at arrival time in each host
          current_community<-data.frame(rumen_otus[rumen_otus$subject==subject & rumen_otus$t==arrival_time,5:2547])
          current_comm_analysis[row,1:3]<-c(subject,arrival_time,persistence)
          current_comm_analysis[row,4:2547]<-current_community
          }
      
      #remove focal OTU from community
      current_comm_analysis$persistence<-as.numeric(current_comm_analysis$persistence)
      focal<-match(OTU,colnames(current_comm_analysis))
      current_comm_analysis<-current_comm_analysis[,-c(focal)]
  
      #record results of adonis test
      dist<-vegdist(current_comm_analysis[,3:2546], method="bray")
      a<-adonis(dist~persistence,permutations=1000,data=current_comm_analysis)
      persistence_F<-a$aov.tab[1,4]
      persistence_pvalue<-a$aov.tab[1,6]
      
      ##dependence on microbiome composition before arrival (t=-1)
      prior_comm_analysis<-data.frame(matrix(nrow=0,ncol=2547))
      colnames(prior_comm_analysis)=c("subject","arrival_time","persistence",colnames(rumen_otus[4:2547]))
      for (row in seq(1,nrow(data_per_OTU))){
          subject<-as.character(data_per_OTU[row,"subjectID"])
          arrival_time<-data_per_OTU[row,"arrival_time"]
          persistence<-data_per_OTU[row,"persistence"]
        
          #in each host, find the last timepoint taken before arrival & extract community
          comm_subset<-rumen_otus[rumen_otus$subject==subject,]
          comm_subset<-comm_subset[order(comm_subset$t),]
          if (match(arrival_time,comm_subset$t)>1){
              prior_timepoint<-as.numeric(comm_subset[match(arrival_time,comm_subset$t)-1,"t"])
              prior_community<-data.frame(comm_subset[comm_subset$t==prior_timepoint,5:2547])
              prior_comm_analysis[row,1:3]<-c(subject,arrival_time,persistence)
              prior_comm_analysis[row,4:2547]<-prior_community
          }
       }
      
      #remove focal OTU from community
      prior_comm_analysis$persistence<-as.numeric(prior_comm_analysis$persistence)
      focal<-match(OTU,colnames(prior_comm_analysis))
      prior_comm_analysis<-prior_comm_analysis[,-c(focal)]
      prior_comm_analysis<-prior_comm_analysis[!is.na(prior_comm_analysis$arrival_time),]
      
      #record results of adonis test
      if (nrow(prior_comm_analysis[prior_comm_analysis$persistence!=1 & prior_comm_analysis$persistence!=0,])>2){
          dist<-vegdist(prior_comm_analysis[,3:2546], method="bray")
          a<-adonis(dist~persistence,permutations=1000,data=prior_comm_analysis)
          persistence_F_prior<-a$aov.tab[1,4]
          persistence_pvalue_prior<-a$aov.tab[1,6]
        }
      else {persistence_F_prior<-NA; persistence_pvalue_prior<-NA}
    
      #iteratively update data frames
      summary_rumen_OTUs<-rbind(summary_rumen_OTUs,data.frame(OTU,mean_arrival,mean_persistence,sd_persistence,num_hosts,taxonomy,persistence_F,persistence_pvalue,persistence_F_prior,persistence_pvalue_prior))
      rumen_list[[i]]<-data_per_OTU
      i=i+1
    }
  }
                                 
                                 
##DESeq analysis of "partner" OTUs whose abundance at (t=0) or before (t=-1) arrival predicts persistence of "focal" OTU

#to be detectably sensitive to microbiome composition, OTUs need to vary a lot in persistence
#here I am considering the set of OTUs with variation in persistence (SD) above the median level-- but this could be adjusted
                                 
sd_persist_med<-summary(summary_rumen_OTUs$sd_persistence)["Median"]

tmp<-summary_rumen_OTUs[summary_rumen_OTUs$sd_persistence>=sd_persist_med,]
tmp$persistence_padj<-p.adjust(tmp$persistence_pvalue,method="BH")
rumen_sensitive_OTUs_current<-data.frame(tmp[tmp$persistence_padj<0.1,"OTU"])$OTU

tmp$persistence_padj_prior<-p.adjust(tmp$persistence_pvalue_prior,method="BH")
rumen_sensitive_OTUs_prior<-data.frame(tmp[tmp$persistence_padj_prior<0.1,"OTU"])$OTU

##identify partner OTUs within the community at arrival (t=0)
rumen_current_deseq<-data.frame(matrix(nrow=0,ncol=22))

#set up dataframe of community at arrival again (for the small number of sensitive OTUs)
for (focal_OTU in rumen_sensitive_OTUs_current){
  data_per_OTU<-rumen_list[[match(focal_OTU,summary_rumen_OTUs$OTU)]]
  current_comm_analysis<-data.frame(matrix(nrow=0,ncol=2547))
  colnames(current_comm_analysis)=c("subject","arrival_time","persistence",colnames(rumen_otus[4:2547]))
  data_per_OTU<-data_per_OTU[!is.na(data_per_OTU$arrival_time) & !is.na(data_per_OTU$persistence),]
  
  for (row in seq(1,nrow(data_per_OTU))){
    subject<-as.character(data_per_OTU[row,"subjectID"])
    arrival_time<-data_per_OTU[row,"arrival_time"]
    persistence<-data_per_OTU[row,"persistence"]
    current_community<-data.frame(rumen_otus[rumen_otus$subject==subject & rumen_otus$t==arrival_time,5:2547])
    current_comm_analysis[row,1:3]<-c(subject,arrival_time,persistence)
    current_comm_analysis[row,4:2547]<-current_community
    }

  focal<-match(OTU,colnames(current_comm_analysis))
  current_comm_analysis<-current_comm_analysis[,-c(focal)]
  current_comm_analysis$persistence<-as.numeric(current_comm_analysis$persistence)
  
  #convert to phyloseq object
  #1: OTU table (convert to numeric, add 1 to everything so DEseq can run)
  phyloseq_otu_table<-t(current_comm_analysis[,4:2546])
  colnames(phyloseq_otu_table)<-current_comm_analysis$subject
  phyloseq_otu_table<-mutate_all(data.frame(phyloseq_otu_table), function(x) as.numeric(as.character(x)))
  phyloseq_otu_table<-mutate_all(data.frame(phyloseq_otu_table), function(x){x+1})
  phyloseq_otu_table<-as.matrix(phyloseq_otu_table)
  phyloseq_otu_table = otu_table(phyloseq_otu_table, taxa_are_rows = TRUE)
  
  #2: taxonomy (remove OTUs that are not present in the OTU table)
  phyloseq_taxonomy<-data.frame(rumen_taxonomy)
  phyloseq_taxonomy<-phyloseq_taxonomy[!phyloseq_taxonomy$X1%in%setdiff(phyloseq_taxonomy$X1,colnames(current_comm_analysis[,4:2546])),]
  phyloseq_taxonomy<-phyloseq_taxonomy[order(phyloseq_taxonomy$X1),]
  colnames(phyloseq_taxonomy)[2:8]=c("Domain","Phylum","Class","Order","Family","Genus","Species")
  phyloseq_taxonomy$Domain<-substr(phyloseq_taxonomy$Domain,4,nchar(phyloseq_taxonomy$Domain))
  phyloseq_taxonomy$Phylum<-substr(phyloseq_taxonomy$Phylum,4,nchar(phyloseq_taxonomy$Phylum))
  phyloseq_taxonomy$Class<-substr(phyloseq_taxonomy$Class,4,nchar(phyloseq_taxonomy$Class))
  phyloseq_taxonomy$Order<-substr(phyloseq_taxonomy$Order,4,nchar(phyloseq_taxonomy$Order))
  phyloseq_taxonomy$Family<-substr(phyloseq_taxonomy$Family,4,nchar(phyloseq_taxonomy$Family))
  phyloseq_taxonomy$Genus<-substr(phyloseq_taxonomy$Genus,4,nchar(phyloseq_taxonomy$Genus))
  phyloseq_taxonomy$Species<-substr(phyloseq_taxonomy$Species,4,nchar(phyloseq_taxonomy$Species))
  rownames(phyloseq_taxonomy)<-taxa_names(phyloseq_otu_table)
  phyloseq_taxonomy<-tax_table(as.matrix(phyloseq_taxonomy))

  #3: sample info (add rownames)
  phyloseq_samples<-data.frame(current_comm_analysis[,c(1:3)])
  phyloseq_samples$subject<-paste("X",phyloseq_samples$subject,sep="")
  rownames(phyloseq_samples)<-phyloseq_samples[,1]
  phyloseq_samples<-phyloseq_samples[,-c(1)]
  phyloseq_samples<-sample_data(phyloseq_samples)
                                                            
  #create a phyloseq object
  phyloseq_obj<-phyloseq(phyloseq_otu_table,phyloseq_taxonomy,phyloseq_samples)
  deseq_test = phyloseq_to_deseq2(phyloseq_obj, ~persistence)
  deseq_output<-DESeq(deseq_test, test="Wald", fitType="parametric")
  #extract significantly enriched/depleted taxa
  res = results(deseq_output, cooksCutoff = FALSE)
  sigtab = res[which(res$padj < 0.05), ]
  sigtab = cbind(as(sigtab, "data.frame"), as(tax_table(phyloseq_obj)[rownames(sigtab), ], "matrix"))
  sigtab<-sigtab[order(-sigtab$log2FoldChange),]
  sigtab$x<-seq(1,nrow(sigtab))
  
  #add focal OTU info and bind to master dataframe
  taxonomy<-rumen_taxonomy[rumen_taxonomy$X1==focal_OTU,3:8]
  sigtab[,"focal_OTU"]<-focal_OTU
  sigtab[,"focal_Phylum"]<-substr(taxonomy[1],4,nchar(taxonomy[1]))
  sigtab[,"focal_Class"]<-substr(taxonomy[2],4,nchar(taxonomy[2]))
  sigtab[,"focal_Order"]<-substr(taxonomy[3],4,nchar(taxonomy[3]))
  sigtab[,"focal_Family"]<-substr(taxonomy[4],4,nchar(taxonomy[4]))
  sigtab[,"focal_Genus"]<-substr(taxonomy[5],4,nchar(taxonomy[5]))
  sigtab[,"focal_Species"]<-substr(taxonomy[6],4,nchar(taxonomy[6]))
  rumen_current_deseq<-rbind(rumen_current_deseq,sigtab)
}       
        
                                 
##identify partner OTUs within the community before arrival (t=-1)                               
rumen_prior_deseq<-data.frame(matrix(nrow=0,ncol=22))
                                 
#set up dataframe of community at arrival again (for the small number of sensitive OTUs)
for (focal_OTU in rumen_sensitive_OTUs_prior){
    data_per_OTU<-rumen_list[[match(focal_OTU,summary_rumen_OTUs$OTU)]]
    prior_comm_analysis<-data.frame(matrix(nrow=0,ncol=2547))
    colnames(prior_comm_analysis)=c("subject","arrival_time","persistence",colnames(rumen_otus[4:2547]))
    data_per_OTU<-data_per_OTU[!is.na(data_per_OTU$arrival_time) & !is.na(data_per_OTU$persistence),]
  
    for (row in seq(1,nrow(data_per_OTU))){
        subject<-as.character(data_per_OTU[row,"subjectID"])
        arrival_time<-data_per_OTU[row,"arrival_time"]
        persistence<-data_per_OTU[row,"persistence"]
        
        #in each host, find the last timepoint taken before arrival & extract community
        comm_subset<-rumen_otus[rumen_otus$subject==subject,]
        comm_subset<-comm_subset[order(comm_subset$t),]
        if (match(arrival_time,comm_subset$t)>1){
            prior_timepoint<-as.numeric(comm_subset[match(arrival_time,comm_subset$t)-1,"t"])
            prior_community<-data.frame(comm_subset[comm_subset$t==prior_timepoint,5:2547])
            prior_comm_analysis[row,1:3]<-c(subject,arrival_time,persistence)
            prior_comm_analysis[row,4:2547]<-prior_community
        }
      }
      
    #remove focal OTU from community
    prior_comm_analysis$persistence<-as.numeric(prior_comm_analysis$persistence)
    focal<-match(OTU,colnames(prior_comm_analysis))
    prior_comm_analysis<-prior_comm_analysis[,-c(focal)]
    prior_comm_analysis<-prior_comm_analysis[!is.na(prior_comm_analysis$arrival_time),]
    prior_comm_analysis<-prior_comm_analysis[!is.na(prior_comm_analysis$arrival_time),]
  
  #convert to phyloseq object
  #1: OTU table (convert to numeric, add 1 to everything so DEseq can run)
  phyloseq_otu_table<-t(prior_comm_analysis[,4:2546])
  colnames(phyloseq_otu_table)<-prior_comm_analysis$subject
  phyloseq_otu_table<-mutate_all(data.frame(phyloseq_otu_table), function(x) as.numeric(as.character(x)))
  phyloseq_otu_table<-mutate_all(data.frame(phyloseq_otu_table), function(x){x+1})
  phyloseq_otu_table<-as.matrix(phyloseq_otu_table)
  phyloseq_otu_table = otu_table(phyloseq_otu_table, taxa_are_rows = TRUE)
  
  #2: taxonomy (remove OTUs that are not present in the OTU table)
  phyloseq_taxonomy<-data.frame(rumen_taxonomy)
  phyloseq_taxonomy<-phyloseq_taxonomy[!phyloseq_taxonomy$X1%in%setdiff(phyloseq_taxonomy$X1,colnames(current_comm_analysis[,4:2546])),]
  phyloseq_taxonomy<-phyloseq_taxonomy[order(phyloseq_taxonomy$X1),]
  colnames(phyloseq_taxonomy)[2:8]=c("Domain","Phylum","Class","Order","Family","Genus","Species")
  phyloseq_taxonomy$Domain<-substr(phyloseq_taxonomy$Domain,4,nchar(phyloseq_taxonomy$Domain))
  phyloseq_taxonomy$Phylum<-substr(phyloseq_taxonomy$Phylum,4,nchar(phyloseq_taxonomy$Phylum))
  phyloseq_taxonomy$Class<-substr(phyloseq_taxonomy$Class,4,nchar(phyloseq_taxonomy$Class))
  phyloseq_taxonomy$Order<-substr(phyloseq_taxonomy$Order,4,nchar(phyloseq_taxonomy$Order))
  phyloseq_taxonomy$Family<-substr(phyloseq_taxonomy$Family,4,nchar(phyloseq_taxonomy$Family))
  phyloseq_taxonomy$Genus<-substr(phyloseq_taxonomy$Genus,4,nchar(phyloseq_taxonomy$Genus))
  phyloseq_taxonomy$Species<-substr(phyloseq_taxonomy$Species,4,nchar(phyloseq_taxonomy$Species))
  rownames(phyloseq_taxonomy)<-taxa_names(phyloseq_otu_table)
  phyloseq_taxonomy<-tax_table(as.matrix(phyloseq_taxonomy))

  #3: sample info (add rownames)
  phyloseq_samples<-data.frame(current_comm_analysis[,c(1:3)])
  phyloseq_samples$subject<-paste("X",phyloseq_samples$subject,sep="")
  rownames(phyloseq_samples)<-phyloseq_samples[,1]
  phyloseq_samples<-phyloseq_samples[,-c(1)]
  phyloseq_samples<-sample_data(phyloseq_samples)

  #create a phyloseq object
  phyloseq_obj<-phyloseq(phyloseq_otu_table,phyloseq_taxonomy,phyloseq_samples)
  deseq_test = phyloseq_to_deseq2(phyloseq_obj, ~persistence)
  deseq_output<-DESeq(deseq_test, test="Wald", fitType="parametric")
  #extract significantly enriched/depleted taxa
  res = results(deseq_output, cooksCutoff = FALSE)
  sigtab = res[which(res$padj < 0.05), ]
  sigtab = cbind(as(sigtab, "data.frame"), as(tax_table(phyloseq_obj)[rownames(sigtab), ], "matrix"))
  sigtab<-sigtab[order(-sigtab$log2FoldChange),]
  sigtab$x<-seq(1,nrow(sigtab))
  
  #add focal OTU info and bind to master dataframe
  taxonomy<-rumen_taxonomy[rumen_taxonomy$X1==focal_OTU,3:8]
  sigtab[,"focal_OTU"]<-focal_OTU
  sigtab[,"focal_Phylum"]<-substr(taxonomy[1],4,nchar(taxonomy[1]))
  sigtab[,"focal_Class"]<-substr(taxonomy[2],4,nchar(taxonomy[2]))
  sigtab[,"focal_Order"]<-substr(taxonomy[3],4,nchar(taxonomy[3]))
  sigtab[,"focal_Family"]<-substr(taxonomy[4],4,nchar(taxonomy[4]))
  sigtab[,"focal_Genus"]<-substr(taxonomy[5],4,nchar(taxonomy[5]))
  sigtab[,"focal_Species"]<-substr(taxonomy[6],4,nchar(taxonomy[6]))
  rumen_prior_deseq<-rbind(rumen_prior_deseq,sigtab)
}
                                 
##Compare relatedness of DESeq-identified pairs to 1000 iterations of random pairs in the cow rumen microbiome
#How many DESeq-identified pairs come from the same family? (t=-1)
rumen_perm_current<-c()
all_pairs<-combn(as.character(summary_rumen_OTUs$OTU),2) #all possible pairs of OTUs that passed our initial filtering criteria
for (i in seq(1,1000)){
  print(i)
  null_relatedness_data<-data.frame(matrix(nrow=nrow(rumen_current_deseq),ncol=5))
  colnames(null_relatedness_data)<-c("OTU_A","Family_A","OTU_B","Family_B","relatedness")
  sampled_pairs<-all_pairs[,sample(seq(1,ncol(all_pairs)),nrow(rumen_current_deseq))] #randomly sample the same number of pairs as in the observe data

  for (j in seq(1,nrow(rumen_current_deseq))){
    OTU_A<-sampled_pairs[1,j]
    family_A<-as.character(rumen_taxonomy[rumen_taxonomy$X1==OTU_A,6])
    OTU_B<-sampled_pairs[2,j]
    family_B<-as.character(rumen_taxonomy[rumen_taxonomy$X1==OTU_B,6])

    if (!is.na(family_A) & !is.na(family_B) & substr(family_A,4,5)!="" & substr(family_B,4,5)!=""){
    relatedness<-as.numeric(family_A==family_B)
    null_relatedness_data[j,]<-c(OTU_A,family_A,OTU_B,family_B,relatedness)
    }
  }
  null_relatedness_data$relatedness<-as.numeric(null_relatedness_data$relatedness)
  rumen_perm_current<-c(rumen_perm_current,mean(null_relatedness_data$relatedness,na.rm=T))
}

obs_means_neg<-nrow(rumen_current_deseq[rumen_current_deseq$log2FoldChange<0 & !is.na(rumen_current_deseq$focal_Family) & !is.na(rumen_current_deseq$Family) & substr(rumen_current_deseq$focal_Family,4,5)!="" & substr(rumen_current_deseq$Family,4,5)!="" & rumen_current_deseq$focal_Family==rumen_current_deseq$Family,])/nrow(rumen_current_deseq[rumen_current_deseq$log2FoldChange<0 & !is.na(rumen_current_deseq$focal_Family) & !is.na(rumen_current_deseq$Family) & substr(rumen_current_deseq$focal_Family,4,5)!="" & substr(rumen_current_deseq$Family,4,5)!="",])
print(length(rumen_perm_current[rumen_perm_current<=obs_means_neg])/1000)
print(length(rumen_perm_current[rumen_perm_current>=obs_means_neg])/1000)

                                 
obs_means_pos<-nrow(rumen_current_deseq[rumen_current_deseq$log2FoldChange>0 & !is.na(rumen_current_deseq$focal_Family) & !is.na(rumen_current_deseq$Family) & substr(rumen_current_deseq$focal_Family,4,5)!="" & substr(rumen_current_deseq$Family,4,5)!="" & rumen_current_deseq$focal_Family==rumen_current_deseq$Family,])/nrow(rumen_current_deseq[rumen_current_deseq$log2FoldChange>0 & !is.na(rumen_current_deseq$focal_Family) & !is.na(rumen_current_deseq$Family) & substr(rumen_current_deseq$focal_Family,4,5)!="" & substr(rumen_current_deseq$Family,4,5)!="",])
print(length(rumen_perm_current[rumen_perm_current<=obs_means_neg])/1000)
print(length(rumen_perm_current[rumen_perm_current>=obs_means_neg])/1000)

                                 
#How many DESeq-identified pairs come from the same family? (t=-1)
rumen_perm_prior<-c()
all_pairs<-combn(as.character(summary_rumen_OTUs$OTU),2) #all possible pairs of OTUs that passed our initial filtering criteria
for (i in seq(1,1000)){
  print(i)
  null_relatedness_data<-data.frame(matrix(nrow=nrow(rumen_prior_deseq),ncol=5))
  colnames(null_relatedness_data)<-c("OTU_A","Family_A","OTU_B","Family_B","relatedness")
  sampled_pairs<-all_pairs[,sample(seq(1,ncol(all_pairs)),nrow(rumen_prior_deseq))] #randomly sample the same number of pairs as in the observe data

  for (j in seq(1,nrow(rumen_prior_deseq))){
    OTU_A<-sampled_pairs[1,j]
    family_A<-as.character(rumen_taxonomy[rumen_taxonomy$X1==OTU_A,6])
    OTU_B<-sampled_pairs[2,j]
    family_B<-as.character(rumen_taxonomy[rumen_taxonomy$X1==OTU_B,6])

    if (!is.na(family_A) & !is.na(family_B) & substr(family_A,4,5)!="" & substr(family_B,4,5)!=""){
    relatedness<-as.numeric(family_A==family_B)
    null_relatedness_data[j,]<-c(OTU_A,family_A,OTU_B,family_B,relatedness)
    }
  }
  null_relatedness_data$relatedness<-as.numeric(null_relatedness_data$relatedness)
  rumen_perm_prior<-c(rumen_perm_prior,mean(null_relatedness_data$relatedness,na.rm=T))
}
                                 
obs_means_neg<-nrow(rumen_prior_deseq[rumen_prior_deseq$log2FoldChange<0 & !is.na(rumen_prior_deseq$focal_Family) & !is.na(rumen_prior_deseq$Family) & substr(rumen_prior_deseq$focal_Family,4,5)!="" & substr(rumen_prior_deseq$Family,4,5)!="" & rumen_prior_deseq$focal_Family==rumen_prior_deseq$Family,])/nrow(rumen_prior_deseq[rumen_prior_deseq$log2FoldChange<0 & !is.na(rumen_prior_deseq$focal_Family) & !is.na(rumen_prior_deseq$Family) & substr(rumen_prior_deseq$focal_Family,4,5)!="" & substr(rumen_prior_deseq$Family,4,5)!="",])
print(length(rumen_perm_prior[rumen_perm_prior<=obs_means_neg])/1000)
print(length(rumen_perm_prior[rumen_perm_prior>=obs_means_neg])/1000)

                                 
obs_means_pos<-nrow(rumen_prior_deseq[rumen_prior_deseq$log2FoldChange>0 & !is.na(rumen_prior_deseq$focal_Family) & !is.na(rumen_prior_deseq$Family) & substr(rumen_prior_deseq$focal_Family,4,5)!="" & substr(rumen_prior_deseq$Family,4,5)!="" & rumen_prior_deseq$focal_Family==rumen_prior_deseq$Family,])/nrow(rumen_prior_deseq[rumen_prior_deseq$log2FoldChange>0 & !is.na(rumen_prior_deseq$focal_Family) & !is.na(rumen_prior_deseq$Family) & substr(rumen_prior_deseq$focal_Family,4,5)!="" & substr(rumen_prior_deseq$Family,4,5)!="",])
print(length(rumen_perm_prior[rumen_perm_prior<=obs_means_pos])/1000)
print(length(rumen_perm_prior[rumen_perm_prior>=obs_means_pos])/1000)

