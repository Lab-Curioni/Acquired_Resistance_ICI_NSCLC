# query CDGS, extract the number of cases harbouring a certain mutation
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

#myTable <- read.table('C:/Users/fzickman/Desktop/cBioportalSNV/TB2_combination_details.txt_listForCbioportalQuery.txt',header=TRUE,sep="\t")
#typeOfCancer = 'Melanoma'
#outfile = 'C:/Users/fzickman/Desktop/cBioportalSNV/cBioportal_TB2_snv'

myTable <- read.table(inputTable,header=TRUE,sep="\t")

positions <- as.numeric(as.character(myTable$Position))  
mygenes <- as.character(myTable$Gene)

print(length(positions))
print(length(mygenes))
print(paste("Cancer type:",typeOfCancer,sep=" "))

outfileMax = paste(outfile,"maxPerVariant.txt",sep="_")
outfileSpecific = paste(outfile,"cancerTypeSpecific.txt",sep="_")
#outfileAll = paste(outfile,"allQueries.txt",sep="_")

if (length(mygenes) == 0){
	writeString = c("study","allCases","gene","position","casesWithMutatedGene","casesWithVariant","frequencyMutatedGene","frequencyVariant")
	write(writeString,file=outfile, ncolumns=length(writeString), sep='\t')
	write(writeString,file=outfileMax, ncolumns=length(writeString), sep='\t')
	write(writeString,file=outfileSpecific, ncolumns=length(writeString), sep='\t')
	print("No genes for query.")
	quit()
}

# for testing: 
#gene = "BRAF"
#mygenes = c("BRAF","EGFR","MIR6859-2")
#position = 140453136
#positions = c(140453136,55229236,14815)
#typeOfCancer = "Melanoma"

#install.packages('cgdsr','/nfs/nas21.ethz.ch/nas/fs2101/eth_nexus_pht_2/utilities/Rlibs/')
library(cgdsr)
mycgds = CGDS("http://www.cbioportal.org/public-portal/")
#test(mycgds)

# get cancer studies
allstudies = getCancerStudies(mycgds)[,1:2]  # first column: study id, second column: description

print(allstudies[,2])

# get ids of cancer type-specific studies, these are always reported
cancerSpecific_ids = grep(typeOfCancer, allstudies[,2], value = F)   # returns the id
cancerSpecific_studies = allstudies[cancerSpecific_ids,]
print(cancerSpecific_studies)

# go through all studies and query several mutations, for large lists, split

studyFreqList = NULL
usedStudies = 0
allCancerTypecases_num = 0
for (i in 1:dim(allstudies)[1])  # for (studyID in melStudies[,1])
{
  #print(i)
  print(paste(i,allstudies[i,1],sep=":"))
  #getProfile
  myProfile = getGeneticProfiles(mycgds,allstudies[i,1])
  # find correct profile id
  myProfile_id = subset(myProfile, myProfile$genetic_alteration_type == 'MUTATION_EXTENDED')$genetic_profile_id
  myProfile_index = 1   # necessary to get correct case list, whcih has a different id
  if (length(myProfile_id) > 1)
  {
    print("Warning! Profile not unique! Using first one!")
    print(myProfile_id)
    myProfile_id = myProfile_id[1]
    print(myProfile_id)
  }
  if (length(myProfile_id) == 0)
  {
    print("Warning! No suitable profile found!")
    #print(myProfile)
    next
  }
  #myProfile_index = which(myProfile$genetic_alteration_type == 'MUTATION_EXTENDED')
  #getCaseList based on correct profile id
  #print(allstudies[i,1])
  myCase = getCaseLists(mycgds,allstudies[i,1])
  numAllCases = length(unlist(strsplit(myCase[grep("Sequenced Tumors", myCase$case_list_name, value = F,ignore.case = TRUE),5], " ", fixed = TRUE, perl = FALSE, useBytes = FALSE)))
  #print(numAllCases)
  if(numAllCases == 0)
  {
    print("Warning! Did not find identifier 'Sequenced Tumors'! Use All Tumors instead!")
    numAllCases = length(unlist(strsplit(myCase[grep("All Tumors", myCase$case_list_name, value = F,ignore.case = TRUE),5], " ", fixed = TRUE, perl = FALSE, useBytes = FALSE)))
    print(numAllCases)
    if(numAllCases == 0)
    {
      print("Warning! Did not find this identifier either! Check study!")
      next
    }
  }
  if(allstudies[i,1] %in% t(cancerSpecific_studies[,1]))
  {
	print(allstudies[i,1])
	allCancerTypecases_num = allCancerTypecases_num + numAllCases
	print(allCancerTypecases_num)
  }
  #getMutation
  usedStudies = usedStudies + 1
  
    for (partPos in 1:ceiling(length(mygenes)/80)) { # Process symbols in chunks of 1000
      start <- 1 + (partPos - 1) * 80
      end <- min(start + 79, length(mygenes))
      #print(paste("Query genes",start,"to",end,sep=" "))
	  Sys.sleep(2)
	  hasFailed <- F
      myMutData = tryCatch(getMutationData(mycgds,allstudies[i,1],myProfile_id,genes=mygenes[start:end]),
	  warning = function(war) {
		print(paste("Warning: No entry for study",allstudies[i,1],sep=" "))
		print(paste("MY_WARNING:  ",war))
	  },
	  error = function(err) {
	  print(paste("Error: No entry for study",allstudies[i,1],sep=" "))
	  print(paste("MY_ERROR:  ",err))
	  hasFailed <- T
	  })
      if (hasFailed)
      {
        #print(myMutData)
        print(paste("No entry for study",allstudies[i,1],sep=" "))
        next
      }
      
      for (j in start:end)
      {
        subsetAll = subset(myMutData, myMutData$gene_symbol == mygenes[j])
		temp = data.frame(subsetAll[!duplicated(subsetAll$case_id), ])
		numCasesMutated = length(temp$gene_symbol)
        #numCasesMutated = length(subsetAll[!duplicated(subsetAll$case_id), ]$gene_symbol)
        #numCasesMutated = length(subset(myMutData, myMutData$gene_symbol == mygenes[j])$gene_symbol)
        subsetVariant = subset(myMutData, myMutData$gene_symbol == mygenes[j] & myMutData$start_position == positions[j])
		tempVariant = data.frame(subsetVariant[!duplicated(subsetVariant$case_id), ])
		numCasesWithVariant = length(tempVariant$gene_symbol)
        #numCasesWithVariant = length(subsetVariant[!duplicated(subsetVariant$case_id), ]$gene_symbol)
        #numCasesWithVariant = length(subset(myMutData, myMutData$gene_symbol == mygenes[j] & myMutData$start_position == positions[j])$gene_symbol)
        freqgenes = (numCasesMutated)/(numAllCases)
        freqVariants = (numCasesWithVariant)/(numAllCases)
        #caseListVector = c(allstudies[i,1],numAllCases,mygenes[j],positions[j],numCasesMutated,numCasesWithVariant,freqgenes,freqVariants)
        studyFreqList = rbind(studyFreqList,c(allstudies[i,1],numAllCases,mygenes[j],positions[j],numCasesMutated,numCasesWithVariant,freqgenes,freqVariants))
      }
    }
    
    rm(partPos, start, end)
  
  
}

print(paste("Used",usedStudies,"studies out of",dim(allstudies)[1],"available studies."),sep=" ")
colnames(studyFreqList) <- c("study","allCases","gene","position","casesWithMutatedGene","casesWithVariant","frequencyMutatedGene","frequencyVariant")
#print(studyFreqList)
studyFreqList_df = data.frame(studyFreqList,row.names = c(1:dim(studyFreqList)[1]))
write.table(studyFreqList_df,file=outfile, quote = F, sep = "\t", row.names = F)
#studyFreqList_df[1:3,]
#find the maximum for each gene, report corresponding study
#print(studyFreqList_df)
maxStudies = NULL
print("MAX_STUDIES")
for (j in 1:length(mygenes))
{
  allCases_gene = data.frame(subset(studyFreqList_df, studyFreqList_df$gene == mygenes[j]))
  allCases_position = data.frame(subset(studyFreqList_df, studyFreqList_df$gene == mygenes[j] & studyFreqList_df$position == positions[j]))
  #print(allCases_gene)
  #print(allCases_position)
  rowNumGene = which.max(as.numeric(as.character(allCases_gene$frequencyMutatedGene)))
  #print(rowNumGene)
  #print(max(as.numeric(as.character(allCases_gene$frequencyMutatedGene))))
  #print(allCases_gene[rowNumGene,])
  rowNumVariant = which.max(as.numeric(as.character(allCases_position$frequencyVariant)))
  #print(allCases_position[rowNumVariant,])
  maxStudies = rbind(maxStudies,allCases_gene[rowNumGene,],allCases_position[rowNumVariant,])
}
colnames(maxStudies) <- c("study","allCases","gene","position","casesWithMutatedGene","casesWithVariant","frequencyMutatedGene","frequencyVariant")
#print(maxStudies)
write.table(maxStudies,file=outfileMax, quote = F, sep = "\t", row.names = F)

# add up for cancer type-specific studies
print(paste("Total cases for specific cancer-type:",allCancerTypecases_num,sep=" "))
allCancerTypeCases = studyFreqList_df[studyFreqList_df$study %in% t(cancerSpecific_studies[1]),]
# only sum up for one gene and position to avoid double counting
#allCancerTypeCases_unique = subset(allCancerTypeCases, allCancerTypeCases$gene == mygenes[1] & allCancerTypeCases$position == positions[1])
#allCancerTypecases_num = sum(as.numeric(as.character(allCancerTypeCases_unique$allCases)))

mutFreqCancerTypeSpecific = NULL
for (j in 1:length(mygenes))
{
  allCancerTypeCases_gene = data.frame(subset(allCancerTypeCases, allCancerTypeCases$gene == mygenes[j]))
  #print("BEFORE")
  #print(allCancerTypeCases_gene)
  allCancerTypeCases_gene = allCancerTypeCases_gene[!duplicated(allCancerTypeCases_gene$study), ]
  #print("AFTER")
  #print(allCancerTypeCases_gene)
  numMutGene = sum(as.numeric(as.character(allCancerTypeCases_gene$casesWithMutatedGene)))
  allCancerTypeCases_position = subset(allCancerTypeCases, allCancerTypeCases$gene == mygenes[j] & allCancerTypeCases$position == positions[j])
  numMutPosition = sum(as.numeric(as.character(allCancerTypeCases_position$casesWithVariant)))
  freqGene = (numMutGene)/(allCancerTypecases_num)
  freqVar = (numMutPosition)/(allCancerTypecases_num)
  mutFreqCancerTypeSpecific = rbind(mutFreqCancerTypeSpecific,c(typeOfCancer,allCancerTypecases_num,mygenes[j],positions[j],numMutGene,numMutPosition,freqGene,freqVar))
}
colnames(mutFreqCancerTypeSpecific) <- c("study","allCases","gene","position","casesWithMutatedGene","casesWithVariant","frequencyMutatedGene","frequencyVariant")
#print(mutFreqCancerTypeSpecific)
write.table(mutFreqCancerTypeSpecific,file=outfileSpecific, quote = F, sep = "\t", row.names = F)
