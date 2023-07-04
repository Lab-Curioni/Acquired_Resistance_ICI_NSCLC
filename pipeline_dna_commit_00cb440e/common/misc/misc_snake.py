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

#config = Config(config)

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
                sample = lineSplit[1]
                output.append(sample)
    return output

def getNormalTumorFiles():
    if not 'SAMPLEMAPPING' in globals():
        return ['NOMAPPINGFILE']
    try:
        open(SAMPLEMAPPING, "r")
    except IOError:
        return ['NOMAPPINGFILE']
    output = []
    sampleMap = dict()
    with open(SAMPLEMAPPING, "r") as f:
        for line in f:
            lineSplit = line.strip().split()
            exp = lineSplit[0]
            sample = lineSplit[1]
            sampleType = lineSplit[2]
            tpoint = lineSplit[3]
            if exp not in sampleMap.keys():
                sampleMap[exp] = dict()
            if tpoint not in sampleMap[exp].keys():
                sampleMap[exp][tpoint] = dict()
            if sampleType not in sampleMap[exp][tpoint].keys():
                sampleMap[exp][tpoint][sampleType] = []
            sampleMap[exp][tpoint][sampleType].append(sample)
    for expKey, expValue in sampleMap.items():
        for tpointKey, tpointValue in expValue.items():
            if 'T' in tpointValue and 'N' in tpointValue:
                for sampleTumor in tpointValue['T']:
                    for sampleNormal in tpointValue['N']:
                        output.append(sampleTumor + '_vs_' + sampleNormal)
    return output

def getNormalFiles():
    if not 'SAMPLEMAPPING' in globals():
        return ['NOMAPPINGFILE']
    try:
        open(SAMPLEMAPPING, "r")
    except IOError:
        return ['NOMAPPINGFILE']
    output = []
    sampleMap = dict()
    with open(SAMPLEMAPPING, "r") as f:
        for line in f:
            lineSplit = line.strip().split()
            exp = lineSplit[0]
            sample = lineSplit[1]
            sampleType = lineSplit[2]
            tpoint = lineSplit[3]
            if exp not in sampleMap.keys():
                sampleMap[exp] = dict()
            if tpoint not in sampleMap[exp].keys():
                sampleMap[exp][tpoint] = dict()
            if sampleType not in sampleMap[exp][tpoint].keys():
                sampleMap[exp][tpoint][sampleType] = []
            sampleMap[exp][tpoint][sampleType].append(sample)
    for expKey, expValue in sampleMap.items():
        for tpointKey, tpointValue in expValue.items():
            if 'T' in tpointValue and 'N' in tpointValue:
                for sampleTumor in tpointValue['T']:
                    for sampleNormal in tpointValue['N']:
                        output.append(sampleNormal)
    return output


def replaceResource(line):
    out=''
    splitValues=line.split('$$$')
    for i in range(0, len(splitValues)):
        if not i%2:
         
            out += splitValues[i]
        else:
            location = splitValues[i].split(':')
            if len(location) == 2:
                try:
                    out += config[location[0]][location[1]]
                except KeyError:
                    return ''
            if len(location) == 3:
                try:
                    if location[1] == "general":
                        out += config[location[0]][location[1]][location[2]]
                    else:
                        out += config[location[0]][ORGANISM][location[2]]
                except KeyError:
                    return ''
    return out

def replaceResourceRecursive(dictionary):
    for key, value in dictionary.items():
        if isinstance(value, dict):
            replaceResourceRecursive(value)
        else:
            dictionary[key] = replaceResource(value)

def postProcessConfigMap():
    # complete all file path in the dictionary
    global TOOLSDIR
    global RESOURCEDIR
    if TOOLSDIR[-1] != '/':
        TOOLSDIR = TOOLSDIR + '/'
    if RESOURCEDIR[-1] != '/':
        RESOURCEDIR += '/'
    config['dirs']={"tools": TOOLSDIR, "resource": RESOURCEDIR}
    for organism in config['resources']:
        for resource in config['resources'][organism]:
            if resource in config['resources']['projectSpecific']:
                config['resources'][organism][resource] = config['resources']['projectSpecific'][resource]
            else:
                config['resources'][organism][resource] = config['dirs']['resource'] + organism + '/' + config['resources'][organism][resource]
    for tool, toolSpecifications in config['tools'].items():
        replaceResourceRecursive(toolSpecifications)

    # GATK specific
    # if not all databases are present adapt the config map
    #if not os.path.isfile(config['resources'][ORGANISM]['Mills_indels']):
    #    print("WARNING: ", config['resources'][ORGANISM]['Mills_indels'], " not present. GATK RealignTargetCreator and GATK BaseRecalibration will not be able to use it.")
    #    config['tools']['GATK']['realignTargetCreator']['known1'] = ''
    #    config['tools']['GATK']['baseRecalibrator']['known1'] = ''
    #if not os.path.isfile(config['resources'][ORGANISM]['1000G_indels']):
    #    print("WARNING: ", config['resources'][ORGANISM]['1000G_indels'], " not present. GATK RealignTargetCreator and GATK BaseRecalibration will not be able to use it.")
    #    config['tools']['GATK']['realignTargetCreator']['known2']= ''
    #    config['tools']['GATK']['baseRecalibrator']['known2'] = ''
    #if not os.path.isfile(config['resources'][ORGANISM]['dbSNP']):
    #    print("WARNING: ", config['resources'][ORGANISM]['dbSNP'], " not present. GATK BaseRecalibration and GATK SNPRecalibrationModel will not be able to use it.")
    #    config['tools']['GATK']['baseRecalibrator']['known3'] = ''
    #    config['tools']['GATK']['gatkSNPrecalibrateModel']['resource4'] = ''
    #if not os.path.isfile(config['resources'][ORGANISM]['hapmap']):
    #    print("WARNING: ", config['resources'][ORGANISM]['hapmap'], " not present. GATK SNPRecalibrationModel will not be able to use it.")
    #    config['tools']['GATK']['gatkSNPrecalibrateModel']['resource1'] = ''
    #if not os.path.isfile(config['resources'][ORGANISM]['1000G_omni']):
    #    print("WARNING: ", config['resources'][ORGANISM]['1000G_omni'], " not present. GATK SNPRecalibrationModel will not be able to use it.")
    #    config['tools']['GATK']['gatkSNPrecalibrateModel']['resource2'] = ''
    #if not os.path.isfile(config['resources'][ORGANISM]['1000G_snp']):
    #    print("WARNING: ", config['resources'][ORGANISM]['1000G_snp'], " not present. GATK SNPRecalibrationModel will not be able to use it.")
    #    config['tools']['GATK']['gatkSNPrecalibrateModel']['resource3'] = ''
    #if not os.path.isfile(config['resources']['general']['gatkKey']):
    #    print("WARNING: GATK will upload a report for each tool used to the Amazon cloud because no GATK key was provided!")
    #    config['tools']['GATK']['call'] = config['tools']['GATK']['call'].strip().split(" --gatk_key")[0]

#rule gunzip:
#    input: 
#        '{sample}.gz'
#    output: 
#        '{sample}'
#    params:
#        lsfoutfile = '{sample}.gunzip.lsfout.log',
#        lsferrfile = '{sample}.gunzip.lsferr.log',
#        scratch = config['tools']['gunzip']['scratch'],
#        mem = config['tools']['gunzip']['memory'],
#        time = config['tools']['gunzip']['time']
#    shell:
#        'gunzip {input}'
