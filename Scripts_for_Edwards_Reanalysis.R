#Analysis of temporal data from rice root endosphere microbiomes
#Robin Herbert

#Load the packages used in this script
library(ggplot2)
library(dplyr)
library(vegan)
library(phyloseq)

#Loads OTU count data, sampling metadata, and the taxonomy of the previously identified OTUs of interest.

otu_Counts<-read.delim("/Users/robinherbert/Desktop/Priority_Effects_Review/Edwards_2018/Sequencing_Reanalysis/Counts.tsv")
rownames(otu_Counts)<-otu_Counts[1:30916,1512]
otu_Counts<-otu_Counts[,-1512]
plant_MetaData<-read.table("/Users/robinherbert/Desktop/Priority_Effects_Review/Edwards_2018/Sequencing_Reanalysis/lc_study_mapping_file.tsv", header = TRUE, row.names = 1)
otu_MetaData<-read.csv("/Users/robinherbert/Desktop/Priority_Effects_Review/Edwards_2018/Sequencing_Reanalysis/predictiveOTUs.csv")

#Limit sampling metadata to that relevant to this reanalysis

new_MetaData<-data.frame(matrix(,nrow = 7, ncol = 1511))
new_MetaData[1,1:ncol(plant_MetaData)]<-as.numeric(plant_MetaData[4,1:ncol(plant_MetaData)])
new_MetaData[2,1:ncol(plant_MetaData)]<-as.character(plant_MetaData[5,1:ncol(plant_MetaData)])
new_MetaData[3,1:ncol(plant_MetaData)]<-as.numeric(plant_MetaData[8,1:ncol(plant_MetaData)])
new_MetaData[4,1:ncol(plant_MetaData)]<-as.character(plant_MetaData[9,1:ncol(plant_MetaData)])
new_MetaData[5,1:ncol(plant_MetaData)]<-as.character(plant_MetaData[10,1:ncol(plant_MetaData)])
new_MetaData[6,1:ncol(plant_MetaData)]<-as.character(plant_MetaData[11,1:ncol(plant_MetaData)])
new_MetaData[7,1:ncol(plant_MetaData)]<-as.numeric(plant_MetaData[12,1:ncol(plant_MetaData)])
colnames(new_MetaData)<-colnames(plant_MetaData)
rownames(new_MetaData)<-c("Age","Compartment","Replicate","Cultivar","State","Site","Season")

#Define lists of timepoints relevant to the reanalysis. All biweekly samples in Edwards, et al., 2018, are included as these are common across all field trials described therein.

plantAge<-c(14,28,42,56,70,84,98,112,126)

#This function makes a list of samples corresponding to a specific root fraction, sampled at the defined timepoints, from each of the three field trials described.

getSampleList<-function(Timepoints, Root_Fraction){
  sampleNumber1<-0
  currentSample<-list()
  sampNames<-c()
  
  for (i in 1:ncol(new_MetaData)){
    if ("2014" %in% new_MetaData[7,i]==TRUE){
      for (j in 1:length(Timepoints)){
        if (Timepoints[j] %in% new_MetaData[1,i]==TRUE){
          if (Root_Fraction %in% new_MetaData[2,i]==TRUE){
            sampleNumber1<-sampleNumber1+1
            currentSample<-colnames(new_MetaData[i])
            sampNames[sampleNumber1]<-currentSample
          }
        }
      }
    }
    else if ("2015" %in% new_MetaData[7,i]==TRUE){
      for (j in 1:length(Timepoints)){
        if (Timepoints[j] %in% new_MetaData[1,i]==TRUE){
          if (Root_Fraction %in% new_MetaData[2,i]==TRUE){
            sampleNumber1<-sampleNumber1+1
            currentSample<-colnames(new_MetaData[i])
            sampNames[sampleNumber1]<-currentSample
          }
        }
      }
    }
    else if ("2016" %in% new_MetaData[7,i]==TRUE){
      for (j in 1:length(Timepoints)){
        if (Timepoints[j] %in% new_MetaData[1,i]==TRUE){
          if (Root_Fraction %in% new_MetaData[2,i]==TRUE){
            sampleNumber1<-sampleNumber1+1
            currentSample<-colnames(new_MetaData[i])
            sampNames[sampleNumber1]<-currentSample
          }
        }
      }
    }
  }
  return(sampNames)
}

currentSampleList<-getSampleList(plantAge,"Endosphere")

#This function calculates relative abundances from the count data representing the above samples

getRelativeAbundanceData<-function(Sample_Names){
  currentSample<-list()
  sampNames<-Sample_Names
  otuTable<-data.frame(matrix(nrow = nrow(otu_Counts),ncol = length(sampNames)))
  rownames(otuTable)<-rownames(otu_Counts)
  colnames(otuTable)<-sampNames
  
  for (j in 1:length(sampNames)){
    current_sample<-match(sampNames[j],colnames(otu_Counts))
    current_community<-as.numeric(otu_Counts[1:nrow(otu_Counts),current_sample])
    community_abundances<-data.frame(matrix(,nrow = 30916,ncol = 1))
    rownames(community_abundances)<-rownames(otu_Counts)
    community_abundances[]<-current_community
    totalCounts<-sum(as.numeric(current_community))
    
    for (k in 1:nrow(community_abundances)){
      relAbund<-as.numeric(current_community[k])/totalCounts
      community_abundances[k,1]<-relAbund
    }
    
    otuTable[,j]<-community_abundances
  }
  
  print(otuTable)
}

relativeAbundanceData<-getRelativeAbundanceData(currentSampleList)

#Define OTUs of interest as those previously described as predictive of plant age in Edwards, et al., 2018

EndosphereOTUs<-otu_MetaData[otu_MetaData$Compartment=="Endosphere",2]

#This function calculates Bray-Curtis dissimilarities between samples and runs nested a PERMANOVA to test if the relative abundance of any of a list of candidate OTUs correlates with community composition, accounting for plant age at sampling

getSampleDissimilarity<-function(OTU_List, Sample_Data){
  
  abundanceTable<-as.matrix(Sample_Data)
  rownames(abundanceTable)<-rownames(Sample_Data)
  colnames(abundanceTable)<-colnames(Sample_Data)
  
  for (j in 1:length(OTU_List)){
    current_mData<-data.frame(matrix(,nrow = ncol(Sample_Data),ncol = 8))
    print(OTU_List[j])
    print(" ")
    current_otu<-match(OTU_List[j],rownames(abundanceTable))
    
    for (i in 1:length(currentSampleList)){
      current_sample<-match(currentSampleList[i],colnames(abundanceTable))
      current_mData[i,8]<-abundanceTable[current_otu,current_sample]
      current_sample<-match(currentSampleList[i],colnames(new_MetaData))
      current_mData[i,1:7]<-new_MetaData[1:7,current_sample]
    }
    
    colnames(current_mData)<-c("Plant_Age","Compartment","Replicate", "Cultivar","State","Field_Site","Year","Rel_Abund_OTUi")
    rownames(current_mData)<-as.character(currentSampleList)
    
    OTU<-otu_table(abundanceTable,taxa_are_rows = TRUE)
    mDATA<-sample_data(current_mData)
    phyObject<-phyloseq(OTU,mDATA)
    
    testDist<-distance(phyObject,type="samples", method = "bray")
    p<-adonis(testDist ~ Plant_Age / Rel_Abund_OTUi, data = data.frame(current_mData))
    print(p)
  }
}

sampleDissimilarity<-getsampleDissimilarity(EndosphereOTUs, relativeAbundanceData)

#This function runs PCoA grouping the samples by those with the top 10% relative abundance of a given candidate OTU

getOrdinationData<-function(OTU_List, Sample_Data){
  
  abundanceTable<-as.matrix(Sample_Data)
  currentAbundanceTable<-matrix(,nrow = nrow(abundanceTable), ncol = ncol(abundanceTable))
  currentAbundanceTable[1:nrow(currentAbundanceTable),1:ncol(currentAbundanceTable)]<-as.numeric(abundanceTable)
  rownames(currentAbundanceTable)<-rownames(abundanceTable)
  colnames(currentAbundanceTable)<-colnames(abundanceTable)
  
  Final_mData<-data.frame(matrix(, ncol = 9))
  colnames(Final_mData)<-c("Plant_Age","Compartment","Replicate", "Cultivar","State","Field_Site","Year","Rel_Abund_OTUi", "Abundance_Percentile")
  rownames(abundanceTable)<-rownames(Sample_Data)
  colnames(abundanceTable)<-colnames(Sample_Data)
  percentile<-10
  
  for (j in 1:length(OTU_List)){
    current_mData<-data.frame(matrix(,nrow = ncol(Sample_Data),ncol = 9))
    print(OTU_List[j])
    print(" ")
    current_otu<-match(OTU_List[j],rownames(abundanceTable))
    
    for (i in 1:length(currentSampleList)){
      current_sample<-match(currentSampleList[i],colnames(abundanceTable))
      current_mData[i,8]<-abundanceTable[current_otu,current_sample]
      current_sample<-match(currentSampleList[i],colnames(new_MetaData))
      current_mData[i,1:7]<-new_MetaData[1:7,current_sample]
    }
    
    colnames(current_mData)<-c("Plant_Age","Compartment","Replicate", "Cultivar","State","Field_Site","Year","Rel_Abund_OTUi", "Abundance_Percentile")
    rownames(current_mData)<-as.character(currentSampleList)
    
    current_mData<-current_mData[order(current_mData$Rel_Abund_OTUi, decreasing = TRUE),]
    
    for (j in 1:length(plantAge)){
      percentile_counter<-0
      currentPlantAge<-current_mData[current_mData$Plant_Age == plantAge[j],]
      
      for (k in 1:nrow(currentPlantAge)){
        percentile_counter<-percentile_counter + 1
        
        if (percentile_counter <= (nrow(currentPlantAge)/percentile)){
          currentPlantAge[k,9]<-"Top 10th Percentile"
        }
        if (percentile_counter > (nrow(currentPlantAge)/percentile)){
          currentPlantAge[k,9]<-"Lower 90th Percentile"
        }
      }
      
      Final_mData<-rbind(Final_mData,currentPlantAge)
    }
    
    colnames(Final_mData)<-c("Plant_Age","Compartment","Replicate", "Cultivar","State","Field_Site","Year","Rel_Abund_OTUi", "Abundance_Percentile")
    Final_mData<-Final_mData[-1,]
    
    #Remove the focal OTU from community
    currentAbundanceTample<-abundanceTable[-c(current_otu),]
    
    OTU<-otu_table(currentAbundanceTable,taxa_are_rows = TRUE)
    mDATA<-sample_data(Final_mData)
    phyObject<-phyloseq(OTU,mDATA)
    forOrdination<-ordinate(phyObject, method = "PCoA", distance = "bray")
    
    p<-plot_ordination(phyObject,forOrdination,type="samples", color = "Abundance_Percentile") + geom_point(size = 3) + ggtitle("Clustering_Rhodoferax_Rich_Samples") + scale_color_manual(name = "Relative Abundance",values = c("black","red")) + theme_classic(base_size=18) + xlab("PCoA Axis 1") + ylab("PCoA Axis 2") + ggtitle("OTU4453710") + theme(plot.title=element_text(hjust=0.5))
    print(p)
  }
}

#Run ordination function focusing on OTU 4453710, annotated as Rhodoferax sp.

ordinationData<-getOrdinationData(4453710, relativeAbundanceData)

#Then retrieve, run statistics, and plot linear regressions of Rhodoferax and three Geobacter spp. present in the list of OTUs of interest

candidateList<-c("Rhodoferax", "Geobacter1","Geobacter2","Geobacter3")
candidateIDs<-c(4453710,344965,167822,137893)

getCovariabilityData<-function(OTU_IDs, OTU_Names, Sample_Data){
  
  abundanceTable<-as.matrix(Sample_Data)
  rownames(abundanceTable)<-rownames(Sample_Data)
  colnames(abundanceTable)<-colnames(Sample_Data)
  metaDataColumns<-7+length(candidateList)
  current_mData<-data.frame(matrix(,nrow = ncol(Sample_Data),ncol = metaDataColumns))
  colnames(current_mData)<-c("Plant_Age","Compartment","Replicate", "Cultivar","State","Field_Site","Year",OTU_Names[1:length(OTU_Names)])
  rownames(current_mData)<-as.character(currentSampleList)
  
  for (i in 1:length(currentSampleList)){
    current_sample<-match(currentSampleList[i],colnames(new_MetaData))
    current_mData[i,1:7]<-new_MetaData[1:7,current_sample]
    
    current_sample<-match(currentSampleList[i],colnames(abundanceTable))
    otuNUM<-7
    
    for (j in 1:length(OTU_IDs)){
      otuNUM<-otuNUM+1
      current_otu<-match(as.character(OTU_IDs[j]),rownames(abundanceTable))
      current_mData[i,otuNUM]<-abundanceTable[current_otu,current_sample]
    }
  }
  
  p<-plot(x=current_mData$Rhodoferax, y=current_mData$Geobacter1,main = "Relative Abundance: Rhodoferax ~ Geobacter1", xlim = c(0,0.1), ylim = c(0,0.02), xlab = "Rhodoferax Relative Abundance", ylab =  "Geobacter sp. 1 Relative Abundance") + theme_classic(base_size=18)
  p<-plot(x=current_mData$Rhodoferax, y=current_mData$Geobacter2,main = "Relative Abundance: Rhodoferax ~ Geobacter2", xlim = c(0,0.1), ylim = c(0,0.02), xlab = "Rhodoferax Relative Abundance", ylab =  "Geobacter sp. 2 Relative Abundance") + theme_classic(base_size=18)
  p<-plot(x=current_mData$Rhodoferax, y=current_mData$Geobacter3,main = "Relative Abundance: Rhodoferax ~ Geobacter3", xlim = c(0,0.1), ylim = c(0,0.02), xlab = "Rhodoferax Relative Abundance", ylab =  "Geobacter sp. 3 Relative Abundance") + theme_classic(base_size=18)
}

getCovariabilityData(candidateIDs, candidateList, relativeAbundanceData)








