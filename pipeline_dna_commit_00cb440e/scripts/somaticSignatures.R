library('SomaticSignatures')
library('VariantAnnotation')
library("lsa")

args = commandArgs(TRUE)
input_vcf = args[1]
output_file = args[2]
genomeVersion = args[3]  # currently supported: hg19, hg38
signatures_dir = args[4] # 30 signtures here: http://cancer.sanger.ac.uk/cancergenome/assets/signatures_probabilities.txt
tumorName = args[5]

signatures_file = paste0(signatures_dir, "signatures30.Rdata")
signatures_descriptions = paste0(signatures_dir, "signatures30.descriptions.txt")
load(signatures_file)

signature_info = read.table(signatures_descriptions, header=T, stringsAsFactors=F, sep='\t')
rownames(signature_info) = signature_info$Signature

#data(signatures21)

if (genomeVersion == "hg19"){
    print("Use BSgenome.Hsapiens.UCSC.hg19")
    library("BSgenome.Hsapiens.UCSC.hg19")
} else if (genomeVersion == "hg38"){
    print("Use BSgenome.Hsapiens.UCSC.hg38")
    library("BSgenome.Hsapiens.UCSC.hg38")
} else {
    print(paste("Error! Could not match genome version: ", genomeVersion, sep=""))
}

vr <- readVcfAsVRanges(input_vcf, genomeVersion)

# remove indels
la = sapply(alt(vr), nchar)
lr = sapply(ref(vr), nchar)
indels = which(la != lr) # might not be conservative enough

vr_snvs = vr
if (length(indels) > 0){ 
    print("Exclude indels")
    vr_snvs = vr[-indels,]
} else {
    print("No indels present")
}

bsGenomeVersion = NULL

if (genomeVersion == "hg19"){
    print("Use hg19")
    bsGenomeVersion = BSgenome.Hsapiens.UCSC.hg19
} else if (genomeVersion == "hg38"){
    print("Use hg38")
    bsGenomeVersion = BSgenome.Hsapiens.UCSC.hg38
} else {
    print(paste("Error! Could not match genome version: ", genomeVersion, sep=""))
}

sca_motifs = mutationContext(vr_snvs, bsGenomeVersion)
sca_mm = motifMatrix(sca_motifs, normalize=T)
idx_tumor = grep(tumorName, colnames(sca_mm)) #grep("Tumor", colnames(sca_mm))
tumor_freqs = sca_mm[,idx_tumor]
#sort(cosine(tumor_freqs, signatures21))
tumor_signatures = as.data.frame(sort(cosine(tumor_freqs, signatures30), decreasing=T))
colnames(tumor_signatures)[1] = "Cosine_similarity"
# add descriptions
tumor_signatures_annotated = merge(tumor_signatures, signature_info, by=0)
tumor_signatures_annotated = subset(tumor_signatures_annotated, select=-c(Row.names))
o = order(tumor_signatures_annotated$Cosine_similarity, decreasing=T)
write.table(tumor_signatures_annotated[o, c("Signature", "Cosine_similarity", "Cancer_types", "Proposed_aetiology", "Add_mut_features", "Comments")], file=output_file, quote=c(3,4,5,6), row.names=F)



