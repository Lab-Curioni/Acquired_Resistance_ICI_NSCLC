import os.path
import sys
import inspect
import copy

fail_instantly = False

class Error(object):
    def __init__(self, key, name):
        self.__key = key
        self.__name = name

    def __add__(self, other):
        return self

    def __call__(self, wildcards=None):
        sys.exit(
            """
            ===============================================
            You have not specified '{}' for '{}'
            ===============================================
            """.format(self.__key, self.__name))

    def __getitem__(self, value):
        return Error(key=self.__key, name=self.__name)

class Config(object):
    def __init__(self, kwargs, name='Config'):
        self.__name = name
        self.__members = {}
        for (key, value) in kwargs.items():
            if isinstance(value, dict):
                self.__members[key] = Config(kwargs=value, name=key)
            else:
                self.__members[key] = value
    
    def __getitem__(self, key):
        if key in self.__members:
            return self.__members[key]
        else:
            if fail_instantly:
                sys.exit(
                    """
                    ===============================================
                    You have not specified '{}' for '{}'
                    ===============================================
                    """.format(key, self.__name))
            else:
                return Error(key=key, name=self.__name)

config = Config(config)

def getSampleNames():
    output = [] #[samplename.replace(FASTQDIR,'').replace('/','')for samplename in glob.glob(FASTQDIR + '*/')]
    if output == []:
        if not 'SAMPLEMAPPING' in globals():
            return ['NOMAPPINGFILE']
        try:
            open(SAMPLEMAPPING, "r")
        except IOError:
            return ['NOMAPPINGFILE']
        sampleMap = dict()
        with open(SAMPLEMAPPING, "r") as f:
            for line in f:
                lineSplit = line.strip().split()
                sample = lineSplit[0]
                output.append(sample)
    return output



def checkFilesAgainstSampleNames(files, sampleNames):
    finalFiles = []
    for f in files:
        for name in sampleNames:
            if name + "/" == f[0:len(name+"/")]:
                finalFiles.append(f)

    return finalFiles

def getSingleFastqFiles(SAMPLENAMES):
    files = [file.replace(FASTQDIR, '').replace('.fastq.gz','')for file in glob.glob(FASTQDIR + '*/SINGLEEND/*.fastq.gz')]
    if files == []:
        files = [file.replace(FASTQDIR, '').replace('.fastq','')for file in glob.glob(FASTQDIR + '*/SINGLEEND/*.fastq')]

    return checkFilesAgainstSampleNames(files, SAMPLENAMES)


def getPairedFastqFiles(SAMPLENAMES):
    files = [file.replace(FASTQDIR, '').replace('.fastq.gz','')for file in glob.glob(FASTQDIR + '*/PAIREDEND/*R[12].fastq.gz')]
    if files == []:
        files = [file.replace(FASTQDIR, '').replace('.fastq','')for file in glob.glob(FASTQDIR + '*/PAIREDEND/*R[12].fastq')]
   
    return checkFilesAgainstSampleNames(files, SAMPLENAMES)

def getPairedFastqFilesWithoutR(SAMPLENAMES):
    files = [file.replace(FASTQDIR, '').replace('_R1.fastq.gz','')for file in glob.glob(FASTQDIR + '*/PAIREDEND/*_R1.fastq.gz')]
    if files == []:
        files = [file.replace(FASTQDIR, '').replace('_R1.fastq','')for file in glob.glob(FASTQDIR + '*/PAIREDEND/*_R1.fastq')]

    return checkFilesAgainstSampleNames(files, SAMPLENAMES)

