import sys
from datetime import date

today=str(date.today())

project_dir=sys.argv[1]
output_file=project_dir + "/read_stats." + today + ".txt"

#feature_suffix=".STAR.GRCm38.featureCounts.txt.summary"
bam_stats_suffix=".rseqc.bam_stat.txt"
gene_count_suffix=".htseq_counts.txt"
sample_list=sys.argv[2]

samples = []

with open(sample_list) as fid:
    for line in fid:
        samples.append(line.strip())


with open(output_file, 'w') as fid_out:
    fid_out.write("Sample\tRaw_reads\tMapping_reads\tUniquely_mapping_reads\tCounted_gene_features\n")
    
    for sample in samples:
        qc_file = "{project_dir}/stats/{sample}{suffix}".format(project_dir=project_dir, sample=sample, suffix=bam_stats_suffix)
        gene_features_file = "{project_dir}/gene_counts/{sample}{suffix}".format(project_dir=project_dir, sample=sample, suffix=gene_count_suffix)
    
        
        with open(qc_file) as fid:
            raw_reads = -1
            not_aligned = -1
            unique_alignments = -1

            for line in fid:
                line = line.strip()
                if line.startswith("Total records:"):
                    raw_reads = int(line.split(':')[1].strip())

                if line.startswith("Unmapped reads:"):
                    not_aligned = int(line.split(':')[1].strip())
                    continue

                if line.startswith("mapq >= mapq_cut (unique):"):
                    unique_alignments = int(line.split(':')[1].strip())
                    continue
   
            reads_aligned = raw_reads - not_aligned               

        # count features

        gene_features=0
        with open(gene_features_file) as fid:
            for line in fid:
                if not line.startswith("__"):
                   gene_features += int(line.split()[1].strip()) 

        fid_out.write("{sample}\t{raw}\t{mapping}\t{unique}\t{gene_features}\n".format(sample=sample, 
            raw=raw_reads, mapping=reads_aligned, unique=unique_alignments, gene_features=gene_features))


