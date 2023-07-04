
import sys


before=sys.argv[1]
resistance=sys.argv[2]
out=sys.argv[3]

before_snvs=open(before,'r')
resistance_snvs=open(resistance,'r')
out_file = open(out, 'w')

found_genes=[]

print("The following snvs were found only in resistance:")
for line in resistance_snvs:
        if line.startswith("Chromosome"):
                out_file.write("Chromosome\tStart\tEnd\tFrequency\tOption\n")
                continue
        gene="".join(line.split("\t")[0:3])
        before_snvs=open(before,'r')
        found=False
        for ln in before_snvs:
                if ln.startswith("Chromosome"):
                        continue
                gene_before="".join(ln.split("\t")[0:3])
                if gene_before==gene:
                        found=True
                        break
        before_snvs.close()
        chrom=line.split("\t")[0].replace("chr","hs")
        pos=line.split("\t")[1]
        freq=line.split("\t")[3]
        if found==False and not gene in found_genes:
                found_genes.append(gene)
                out_file.write(chrom+"\t"+pos+"\t"+pos+"\t"+freq+"\tcolor=red\n")
        if found==True and not gene in found_genes:
                found_genes.append(gene)
                out_file.write(chrom+"\t"+pos+"\t"+pos+"\t"+freq+"\tcolor=grey\n")
before_snvs=open(before,'r')
found=False
for ln in before_snvs:
        chrom=ln.split("\t")[0].replace("chr","hs")
        pos=ln.split("\t")[1]
        freq=ln.split("\t")[3]
        if ln.startswith("Chromosome"):
                continue
        gene_before="".join(ln.split("\t")[0:3])
        if not gene_before in found_genes:
                out_file.write(chrom+"\t"+pos+"\t"+pos+"\t"+freq+"\tcolor=green\n")
before_snvs.close()

