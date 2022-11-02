import FWCore.ParameterSet.Config as cms
import FWCore.Utilities.FileUtils as FileUtils
import FWCore.ParameterSet.VarParsing as VarParsing


# PART 1 : PARSE ARGUMENTS

options = VarParsing.VarParsing ('analysis')
options.register('tm',
                -1,
                VarParsing.VarParsing.multiplicity.singleton,
                VarParsing.VarParsing.varType.int,
                "Time slice of pattern files")
options.register('effTM',
                -1,
                VarParsing.VarParsing.multiplicity.singleton,
                VarParsing.VarParsing.varType.int,
                "Effective TM of HGCal (i.e. for Correlator TMUX 6, HGCal looks at 1 in 3 events in input file, rather than 1 in 18, so set effTM to 3)")
options.parseArguments()

inputFiles = []
for filePath in options.inputFiles:
    if filePath.endswith(".root"):
        inputFiles.append(filePath)
    elif filePath.endswith("_cff.py"):
        filePath = filePath.replace("/python/","/")
        filePath = filePath.replace("/", ".")
        inputFilesImport = getattr(__import__(filePath.strip(".py"),fromlist=["readFiles"]),"readFiles")
        inputFiles.extend( inputFilesImport )
    else:
        inputFiles += FileUtils.loadListFromFile(filePath)

# PART 2: SETUP MAIN CMSSW PROCESS 

process = cms.Process("Stage2FileWriter")

process.load('Configuration.Geometry.GeometryExtended2026D88Reco_cff')
process.load('Configuration.Geometry.GeometryExtended2026D88_cff')
process.load('Configuration.StandardSequences.MagneticField_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:upgradePLS3', '')
process.load("FWCore.MessageLogger.MessageLogger_cfi")

process.source = cms.Source("PoolSource", fileNames = cms.untracked.vstring(inputFiles) )
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(options.maxEvents) )

process.load('L1Trigger.L1THGCalUtilities.stage2FileWriter_cff')

process.MessageLogger.cerr.FwkReport.reportEvery = 1
process.Timing = cms.Service("Timing", summaryOnly = cms.untracked.bool(True))

process.load('FWCore.Modules.preScaler_cfi')


# Correlator tests
if options.effTM == 3 :
    if options.tm == 0:
        # Time slice 0
        process.preScaler.prescaleFactor = options.effTM
        process.preScaler.prescaleOffset = 1
        process.Stage2FileWriter.tmIndex = options.tm
    elif options.tm == 6:
        # #Time slice 6
        process.preScaler.prescaleFactor = options.effTM
        process.preScaler.prescaleOffset = 2
        process.Stage2FileWriter.tmIndex = options.tm
    elif options.tm == 12:
        # Time slice 12
        process.preScaler.prescaleFactor = options.effTM
        process.preScaler.prescaleOffset = 0
        process.Stage2FileWriter.tmIndex = options.tm
# GCT tests
elif options.effTM == 2 :
    if options.tm == 0:
        process.preScaler.prescaleFactor = options.effTM
        process.preScaler.prescaleOffset = 1
        process.Stage2FileWriter.tmIndex = options.tm
    if options.tm == 9:
        process.preScaler.prescaleFactor = options.effTM
        process.preScaler.prescaleOffset = 0
        process.Stage2FileWriter.tmIndex = options.tm
else:
    print ("Producing pattern files from all time slices...")
    process.preScaler.prescaleFactor = 1
    process.preScaler.prescaleOffset = 0
    process.Stage2FileWriter.tmIndex = 99

process.p = cms.Path(process.preScaler*process.Stage2FileWriter)
