import os, glob, sys

# This function adapts the config object to include full path information
include: "../common/misc/misc_snake.py"
postProcessConfigMap()

# Check if the uses specified the proper input and output directories
if not 'VCFDIR' in globals():  # panel only, we start from vcf
    print('You have to specify the root directory of the vcf files!')
    sys.exit(1)
if not 'OUTDIR' in globals():
    print('You have to specify the root directory where the results will be generated!')
    sys.exit(1)
if not 'TMPDIR' in globals():
    print('You have to specify the root directory where temporary files will be stored!')
    sys.exit(1)

# This is the default order in which the programs are executed
# If the user specified a different order the user specified version is chosen.

if not 'CONVERTPANELVCFIN' in globals():
    CONVERTPANELVCFIN = VCFDIR
if not 'CONVERTPANELVCFOUT' in globals():
    CONVERTPANELVCFOUT = OUTDIR + 'convertedPanelVCFs/'
if not 'PARSEANNOTATEDVCFIN' in globals():
    PARSEANNOTATEDVCFIN = CONVERTPANELVCFOUT
if not 'DATABASEQUERY' in globals():
    DATABASEQUERY = OUTDIR + 'databaseQuery/'
if not 'DOWNLOADCLINICALTRIALSOUT' in globals():
    DOWNLOADCLINICALTRIALSOUT = OUTDIR + 'clinicalTrials/'


# Include the rules
include: "clinical_rules_snake.py"

# This rule defines which files should be created
rule clinical:
    input:
	    expand(DATABASEQUERY + '{tumorNormalMatching}.converted.snpEff.dbSNP.clinVar.cosmic.dbnsfp.overview.txt', tumorNormalMatching = getNormalTumorFiles()),
	    expand(DATABASEQUERY + '{tumorNormalMatching}.converted.snpEff.dbSNP.clinVar.cosmic.dbnsfp.overview.dgidb.txt.CompleteTable.txt', tumorNormalMatching = getNormalTumorFiles()),
	    expand(DATABASEQUERY + '{tumorNormalMatching}.converted.snpEff.dbSNP.clinVar.cosmic.dbnsfp.overview.dgidb.txt.CompleteTable.ClinicalTrials.txt', tumorNormalMatching = getNormalTumorFiles()),
	    expand(DATABASEQUERY + '{tumorNormalMatching}.converted.snpEff.dbSNP.clinVar.cosmic.dbnsfp.overview.databaseQueriesSNV.txt_dgidbIndependent.txt_damaging.txt', tumorNormalMatching = getNormalTumorFiles())
        #expand(CONVERTPANELVCFOUT + '{tumorNormalMatching}.converted.vcf', tumorNormalMatching = getNormalTumorFiles())
    output:
        OUTDIR + 'complete.txt'
    params:
        lsfoutfile = OUTDIR + 'complete.lsfout.log',
        lsferrfile = OUTDIR + 'complete.lsferr.log'
    shell:
        'touch {output}'
