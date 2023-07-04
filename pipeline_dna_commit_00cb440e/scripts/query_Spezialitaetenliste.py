# Script scans drugs in a mtbz .txt file against the Preparations.xml file downloaded from http://www.spezialitaetenliste.ch/. 
# If a drug is contained in the Preparations.xml file and contains the tag Oncologica a description is printed out.
# Also considers entries as valid if they differ in their final letter. (e.g. decitabin for decitabine)

import sys
import xml.etree.ElementTree


if len(sys.argv) < 3:
        print "Search drugs (Oncologica) approved in switzerland."
        print "Usage: python2.7 search_substances_swissmedic.py [vcfFile] [DB] [out]"
        sys.exit(1)

vcfFile = sys.argv[1]
swissmedic = sys.argv[2]
out = sys.argv[3]

drug_db=open(swissmedic,'r')
root = xml.etree.ElementTree.parse(drug_db).getroot()
infile = open(vcfFile,'r')
outfile = open(out, 'w')

for line in infile:
        if line.startswith("Gene"):
                outfile.write(line.strip("\n") + "\tApproved/Name\tType\tDescription\n")
        else:
	    print(line.split("\t")[0])
            drugs=line.split("\t")[1].split(",")
            not_approved=[]
            for Drug in drugs:
                drug=Drug.split(" (")[0].lower()
                zugelassen=False
                description=[]
                for Preparation in root.iter('Preparation'):
                    for DescriptionLa in Preparation.iter('DescriptionLa'):
                        if zugelassen==True:
                            break
                        if DescriptionLa.text.lower().startswith(drug) or DescriptionLa.text.lower().startswith(drug[:-1]):
                            description=[]
                            text=""
                            zugelassen=True
                            name= DescriptionLa.text.lower()
                            IT=Preparation.find('ItCodes')
                            for Description in IT.iter('DescriptionDe'):
                                description.append(Description.text)
                                if Description.text=="Oncologica":
                                    for Limitation in Preparation.iter('Limitation'):
                                        text=Limitation.find('DescriptionDe').text.replace('\n','').rstrip()
                            join_description=";".join(description)
                            outfile.write(line.strip("\n").split("\t")[0] + "\t" + Drug.strip("\n") + "\t" + name + "\t" + join_description.encode('utf-8') + "\t" + text.encode('utf-8') + "\n")
                            break
                        else:
                            continue
                if zugelassen is False:
                    outfile.write(line.strip("\n").split("\t")[0] + "\t" + Drug.strip("\n") + "\tnot_approved\n")
