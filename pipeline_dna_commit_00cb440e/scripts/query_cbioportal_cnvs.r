# query CDGS, extract the number of cases harbouring a certain copy number alteration
# Franziska Singer, March 2016

### Read arguments from command line
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
  stop("Usage: Rscript query_cbioportal.r inputTable typeOfCancer outfile")
} else {
  inputTable <- args[1]
  typeOfCancer <- args[2]
  outfile <- args[3]
}
#typeOfCancer = 'Melanoma'
#outfile = '//nas21.ethz.ch/eth_nexus_pht_2/mtbz/test/includeDetailedInfos/smallCNVtest_RSTUDIO'
#inputTable = '//nas21.ethz.ch/eth_nexus_pht_2/mtbz/test/includeDetailedInfos/input/smallCNVtest.txt'

myTable <- read.table(inputTable,header=TRUE,sep="\t")


mycnvs <- as.numeric(as.character(myTable$Copynumber))  
mygenes <- as.character(myTable$Gene)
print(paste("Cancer type:",typeOfCancer,sep=" "))

outfileMax = paste(outfile,"maxPerVariant.txt",sep="_")
outfileSpecific = paste(outfile,"cancerTypeSpecific.txt",sep="_")
#outfileAll = paste(outfile,"allQueriesCNV.txt",sep="_")

if (length(mygenes) == 0){
	writeString = c("study","allCases","gene","CNA_type","casesWithMutatedGene","casesWithVariant","frequencyMutatedGene","frequencyVariant")
	write(writeString,file=outfile, ncolumns=length(writeString), sep='\t')
	write(writeString,file=outfileMax, ncolumns=length(writeString), sep='\t')
	write(writeString,file=outfileSpecific, ncolumns=length(writeString), sep='\t')
	print("No genes for query.")
	quit()
}


# for testing: 
#gene = "BRAF"
#mygenes = c("BRAF","EGFR","MIR6859-2","GNAQ")
#position = 140453136
#mycnvs = c(2.7,1.5,0.5,3.7)

cnvTypeFunction <- function(cnv){
  if (cnv < 2)
  {
    if (cnv <= 0.5)
    {
      return(-2)
    }else
    {
      return(-1)
    }
  }
  if (cnv == 2)
  {
    return(0)
  }
  if (cnv > 2)
  {
    if (cnv >= 3.5)
    {
      return(2)
    }else
    {
      return(1)
    }
  }
}

# different classes in TCGA data: #Discrete copy-number calls from the RAE algorithm. -2 = putative homozygous deletion; -1 = putative hemizygous deletion; 0 = copy-neutral; 1 = low-level gain; 2 = high-level amplification.
mycnvTypes = t(lapply(mycnvs, cnvTypeFunction))

library(cgdsr)
mycgds = CGDS("http://www.cbioportal.org/public-portal/")
#test(mycgds)

# get cancer studies
allstudies = getCancerStudies(mycgds)[,1:2]  # first column: study id, second column: description

# get ids of cancer type-specific studies, these are always reported
cancerSpecific_ids = grep(typeOfCancer, allstudies[,2], value = F)   # returns the id
cancerSpecific_studies = allstudies[cancerSpecific_ids,]
print(cancerSpecific_studies[,1])
print(cancerSpecific_studies)
# go through all studies and query several mutations

studyFreqListCNV = NULL
usedStudies = 0
allCancerTypecases_num = 0
for (i in 1:dim(allstudies)[1])  # for (studyID in melStudies[,1])
{
  print(paste(i,allstudies[i,1],sep=":"))
  #getCaseList
  myProfile = getGeneticProfiles(mycgds,allstudies[i,1])
  
  myProfile_id = subset(myProfile, myProfile$genetic_alteration_type == 'COPY_NUMBER_ALTERATION' & myProfile$show_profile_in_analysis_tab == 'true')$genetic_profile_id
  
  if (length(myProfile_id) > 1)
  {
    print("Warning! Profile not unique! Searching for 'GISTIC' to be concordant with othe studies!")
    #print(paste(myProfile$genetic_profile_id,myProfile$genetic_alteration_type,myProfile$show_profile_in_analysis_tab,sep=" | "))
    myIDs = grep('GISTIC', myProfile$genetic_profile_description, value = F,ignore.case = TRUE)   # returns the id
    potentialProfiles = myProfile[myIDs,]
    myProfile_id = subset(potentialProfiles, potentialProfiles$genetic_alteration_type == 'COPY_NUMBER_ALTERATION' & potentialProfiles$show_profile_in_analysis_tab == 'true')$genetic_profile_id
    #print(myProfile_id)
    #myProfile_id = myProfile_id[1]
    
  }
  if (length(myProfile_id) == 0)
  {
    #print(paste(myProfile$genetic_profile_id,myProfile$genetic_alteration_type,myProfile$show_profile_in_analysis_tab,sep=" | "))
    print("Warning! No suitable profile found!")
    #print(myProfile)
    next
  }
  
  #getCaseList
  myCase = getCaseLists(mycgds,allstudies[i,1])
  numAllCases = length(unlist(strsplit(myCase[grep("Tumor Samples with CNA data", myCase$case_list_name, value = F,ignore.case = TRUE),5], " ", fixed = TRUE, perl = FALSE, useBytes = FALSE)))
  myCaseList = myCase[grep("Tumor Samples with CNA data", myCase$case_list_name, value = F,ignore.case = TRUE),1]
  #print(numAllCases)
  if(numAllCases == 0)
  {
    print("Warning! Did not find identifier 'Tumor Samples with CNA data'! Use 'TUMORS CNA' instead!")
    #print(paste(myCase$case_list_id,myCase$case_list_name,sep=" | "))
    numAllCases = length(unlist(strsplit(myCase[grep("TUMORS CNA", myCase$case_list_name, value = F,ignore.case = TRUE),5], " ", fixed = TRUE, perl = FALSE, useBytes = FALSE)))
    myCaseList = myCase[grep("TUMORS CNA", myCase$case_list_name, value = F,ignore.case = TRUE),5]
    if(numAllCases == 0)
    {
      #print(paste(myCase$case_list_id,myCase$case_list_name,sep=" | "))
      #print(myCase$case_list_name)
      print("Warning! Did not find this identifier either! Check study!")
      next
    }
  }
  #print(numAllCases)
  #print(paste("MY CASES:",myCaseList,sep=" "))
  usedStudies = usedStudies + 1
  #getMutation
  if(allstudies[i,1] %in% t(cancerSpecific_studies[,1]))
  {
	print(paste("Cancer type specific:",allstudies[i,1],sep=""))
	allCancerTypecases_num = allCancerTypecases_num + numAllCases
	#print(allCancerTypecases_num)
  }
  # count samples with any mutation in this gene
  

    for (partPos in 1:ceiling(length(mygenes)/100)) { # Process symbols in chunks of 100
      start <- 1 + (partPos - 1) * 100
      end <- min(start + 99, length(mygenes))
      #print(paste("Query genes",start,"to",end,sep=" "))
      #print("CHUNK Before")
      #print(mygenes[start:end])
	  Sys.sleep(2)
	  hasFailed <- F
      myProfileData = tryCatch(getProfileData(mycgds,genes=mygenes[start:end],myProfile_id,myCaseList),
	  warning = function(war) {
		print(paste("Warning: No entry for study",allstudies[i,1],sep=" "))
		print(paste("MY_WARNING:  ",war))
	  },
	  error = function(err) {
	  print(paste("Error: No entry for study",allstudies[i,1],sep=" "))
	  print(paste("MY_ERROR:  ",err))
	  hasFailed <- T
	  })
	#print("Chunk After")
      #print(myProfileData)
      if (hasFailed)
      {
		print(paste("Abort query for study",allstudies[i,1],sep=" "))
        next
      }
	  noEntryGene = 0
      for (j in start:end)
      {
        colID = which(colnames(myProfileData) == mygenes[j])
        if(length(colID) == 0)
        {
          #print(paste("No entry for",mygenes[j],sep=" "))
		  noEntryGene = noEntryGene + 1
          next
        }
        #subsetAll = subset(myMutData, myMutData$gene_symbol == mygenes[j])
		#temp = data.frame(subsetAll[!duplicated(subsetAll$case_id), ])
		#numCasesMutated = length(temp$gene_symbol)
		
		#subsetTemp = data.frame(subset[!duplicated(data.frame(myProfileData)[0]), ])
		#numCasesMutated = dim(subset(subsetTemp, subsetTemp[colID] != 0 & subsetTemp[colID] != 'NaN')[colID])[1]
		#numCasesWithVariant = dim(subset(subsetTemp, subsetTemp[colID] == mycnvTypes[j] & subsetTemp[colID] != 'NaN')[colID])[1]
		numCasesMutated = dim(subset(myProfileData, myProfileData[colID] != 0 & myProfileData[colID] != 'NaN')[colID])[1]
		#print(subset(myProfileData, myProfileData[colID] != 0 & myProfileData[colID] != 'NaN'))
        numCasesWithVariant = dim(subset(myProfileData, myProfileData[colID] == mycnvTypes[j] & myProfileData[colID] != 'NaN')[colID])[1] 
        freqgenes = (numCasesMutated)/(numAllCases)
        freqVariants = (numCasesWithVariant)/(numAllCases)
        studyFreqListCNV = rbind(studyFreqListCNV,c(allstudies[i,1],numAllCases,mygenes[j],mycnvTypes[j],numCasesMutated,numCasesWithVariant,freqgenes,freqVariants))
      }
	  print(paste("No entry for",noEntryGene,"genes.",sep=" "))
    }
    
    rm(partPos, start, end)
  
  
}

print(paste("Used",usedStudies,"studies out of",dim(allstudies)[1],"available studies."),sep=" ")
colnames(studyFreqListCNV) <- c("study","allCases","gene","CNA_type","casesWithMutatedGene","casesWithVariant","frequencyMutatedGene","frequencyVariant")
#print(studyFreqListCNV)
studyFreqList_df = data.frame(studyFreqListCNV,row.names = c(1:dim(studyFreqListCNV)[1]))
allStudies.df <- data.frame(lapply(data.frame(studyFreqList_df), as.character), stringsAsFactors=FALSE)
write.table(allStudies.df,file=outfile, quote = F, sep = "\t", row.names = F)
#find the maximum for each gene, report corresponding study
maxStudies = NULL

for (j in 1:length(mygenes))
{
  allCases_gene = data.frame(subset(studyFreqList_df, as.character(studyFreqList_df$gene) == mygenes[j]))
  allCases_type = data.frame(subset(studyFreqList_df, as.character(studyFreqList_df$gene) == mygenes[j] & as.numeric(as.character(studyFreqList_df$CNA_type)) == mycnvTypes[j]))
  #print(allCases_gene)
  #print(allCases_position)
  rowNumGene = which.max(as.numeric(as.character(allCases_gene$frequencyMutatedGene)))
  #print(max(as.numeric(as.character(allCases_gene$frequencyMutatedGene))))
  #print(allCases_gene[rowNumGene,])
  rowNumVariant = which.max(as.numeric(as.character(allCases_type$frequencyVariant)))
  #print(allCases_position[rowNumVariant,])
  maxStudies = rbind(maxStudies,allCases_gene[rowNumGene,],allCases_type[rowNumVariant,])
}
colnames(maxStudies) <- c("study","allCases","gene","CNA_type","casesWithMutatedGene","casesWithVariant","frequencyMutatedGene","frequencyVariant")

max.df <- data.frame(lapply(data.frame(maxStudies), as.character), stringsAsFactors=FALSE)
print("MAX_STUDIES")
#print(max.df)
write.table(max.df,file=outfileMax, quote = F, sep = "\t", row.names = F)

# add up for cancer type-specific studies

print(paste("Total cases for specific cancer-type:",allCancerTypecases_num,sep=" "))
allCancerTypeCases = studyFreqList_df[studyFreqList_df$study %in% t(cancerSpecific_studies[1]),]
#print(allCancerTypeCases)
# only sum up for one gene and position to avoid double counting
#allCancerTypeCases_unique = subset(allCancerTypeCases, as.character(allCancerTypeCases$gene) == as.character(allCancerTypeCases[1,3]) & as.numeric(as.character(allCancerTypeCases$CNA_type)) == as.numeric(as.character(allCancerTypeCases[1,4])))
#allCancerTypecases_num = sum(as.numeric(as.character(allCancerTypeCases_unique$allCases)))

mutFreqCancerTypeSpecific = NULL
for (j in 1:length(mygenes))
{
  allCancerTypeCases_gene = subset(allCancerTypeCases, as.character(allCancerTypeCases$gene) == mygenes[j])
  #print("BEFORE")
  #print(allCancerTypeCases_gene)
  allCancerTypeCases_gene = allCancerTypeCases_gene[!duplicated(allCancerTypeCases_gene$study), ]
  #print("AFTER")
  #print(allCancerTypeCases_gene)
  numMutGene = sum(as.numeric(as.character(allCancerTypeCases_gene$casesWithMutatedGene)))
  allCancerTypeCases_type = subset(allCancerTypeCases, as.character(allCancerTypeCases$gene) == mygenes[j] & as.numeric(as.character(allCancerTypeCases$CNA_type)) == mycnvTypes[j])
  allCancerTypeCases_type = allCancerTypeCases_type[!duplicated(allCancerTypeCases_type$study), ]
  numMutPosition = sum(as.numeric(as.character(allCancerTypeCases_type$casesWithVariant)))
  freqGene = (numMutGene)/(allCancerTypecases_num)
  freqVar = (numMutPosition)/(allCancerTypecases_num)
  mutFreqCancerTypeSpecific = rbind(mutFreqCancerTypeSpecific,c(typeOfCancer,allCancerTypecases_num,mygenes[j],mycnvTypes[j],numMutGene,numMutPosition,freqGene,freqVar))
}
colnames(mutFreqCancerTypeSpecific) <- c("study","allCases","gene","CNA_type","casesWithMutatedGene","casesWithVariant","frequencyMutatedGene","frequencyVariant")
#print(mutFreqCancerTypeSpecific)
cancerType.df <- data.frame(lapply(data.frame(mutFreqCancerTypeSpecific), as.character), stringsAsFactors=FALSE)
write.table(cancerType.df,file=outfileSpecific, quote = F, sep = "\t", row.names = F)
