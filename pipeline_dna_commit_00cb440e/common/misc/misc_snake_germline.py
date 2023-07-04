import os.path

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

def getExperiments():
    output = [] 
    if output == []:
        if not 'SAMPLEMAPPING' in globals():
            return ['NOMAPPINGFILE']
        try:
            open(SAMPLEMAPPING, "r")
        except IOError:
            return ['NOMAPPINGFILE']
        with open(SAMPLEMAPPING, "r") as f:
            for line in f:
                lineSplit = line.strip().split()
                experiment = lineSplit[0]
                output.append(experiment)
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

