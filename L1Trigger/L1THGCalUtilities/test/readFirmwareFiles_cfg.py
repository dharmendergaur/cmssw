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
options.register('sector',
                -1,
                VarParsing.VarParsing.multiplicity.singleton,
                VarParsing.VarParsing.varType.int,
                "Sector of pattern files")
options.parseArguments()

inputFiles = []
for filePath in options.inputFiles:
    if filePath.endswith(".root"):
        inputFiles.append(filePath)
    else:
        inputFiles += FileUtils.loadListFromFile(filePath)


# PART 2: SETUP MAIN CMSSW PROCESS 

process = cms.Process("HGCClusterValidation")

process.load('Configuration.Geometry.GeometryExtended2026D88Reco_cff')
process.load('Configuration.Geometry.GeometryExtended2026D88_cff')
process.load('Configuration.StandardSequences.MagneticField_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:upgradePLS3', '')
process.load("FWCore.MessageLogger.MessageLogger_cfi")

process.source = cms.Source("PoolSource", fileNames = cms.untracked.vstring(inputFiles) )
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(options.maxEvents) )
process.load('L1Trigger.L1THGCalUtilities.stage2FileReader_cff')

process.MessageLogger.cerr.FwkReport.reportEvery = 1
process.Timing = cms.Service("Timing", summaryOnly = cms.untracked.bool(True))

# process.Stage2FileReader.files = cms.vstring("output_hgcAlgo.txt")

process.Stage2FileReader.files = cms.vstring(
    # 'hgc_sec4_tm6-output-hw_0.txt',
    # 'hgc_sec4_tm6-output-hw_1.txt',
    # 'hgc_sec4_tm6-output-hw_2.txt',
    # 'output_hgcAlgo.txt'
)
# for i in range(8): process.Stage2FileReader.files.append(f'hgc_sec4_tm6-output-hw_{i}.txt')
for i in range(95): process.Stage2FileReader.files.append(f'patterns/hgc_sec5_tm99-output-hw-{i}.txt')
process.Stage2FileReader.sector = options.sector


process.load('FWCore.Modules.preScaler_cfi')
# process.preScaler.prescaleFactor = 3
# process.preScaler.prescaleOffset = 0

if options.tm == 0:
    # Time slice 0
    process.preScaler.prescaleFactor = 3
    process.preScaler.prescaleOffset = 1
elif options.tm == 6:
    # #Time slice 6
    process.preScaler.prescaleFactor = 3
    process.preScaler.prescaleOffset = 2
elif options.tm == 12:
    # Time slice 12
    process.preScaler.prescaleFactor = 3
    process.preScaler.prescaleOffset = 0
else:
    print ("Reading pattern files from all time slices...")
    process.preScaler.prescaleFactor = 1
    process.preScaler.prescaleOffset = 0

process.p = cms.Path(process.preScaler*process.Stage2FileReader)

# load ntuplizer
process.load('Configuration.StandardSequences.Services_cff')
process.load('SimGeneral.HepPDTESSource.pythiapdt_cfi')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.EventContent.EventContent_cff')
process.load('SimGeneral.MixingModule.mixNoPU_cfi')
process.load('Configuration.Geometry.GeometryExtended2026D88Reco_cff')
process.load('Configuration.Geometry.GeometryExtended2026D88_cff')
process.load('Configuration.StandardSequences.MagneticField_cff')
process.load('Configuration.StandardSequences.Generator_cff')
process.load('IOMC.EventVertexGenerators.VtxSmearedHLLHC14TeV_cfi')
process.load('GeneratorInterface.Core.genFilterSummary_cff')
process.load('Configuration.StandardSequences.SimIdeal_cff')
process.load('Configuration.StandardSequences.Digi_cff')
process.load('Configuration.StandardSequences.SimL1Emulator_cff')
process.load('Configuration.StandardSequences.DigiToRaw_cff')
process.load('Configuration.StandardSequences.EndOfProcess_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:phase2_realistic_T21', '')

process.TFileService = cms.Service(
    "TFileService",
    fileName = cms.string("ntuple.root")
    )
process.load('L1Trigger.L1THGCalUtilities.hgcalTriggerNtuples_cff')
from L1Trigger.L1THGCalUtilities.customNtuples import custom_ntuples_standalone_clustering, custom_ntuples_standalone_tower
process = custom_ntuples_standalone_clustering(process)
process = custom_ntuples_standalone_tower(process)

process.l1tHGCalTriggerNtuplizer.Ntuples = cms.VPSet( process.ntuple_multiclusters)
# process.ntuple_multiclusters.Multiclusters = ('Stage2FileReader:HWClusters')
process.ntuple_multiclusters.Multiclusters = ('Stage2FileReader:RefClusters')

process.ntuple_step = cms.Path(process.L1THGCalTriggerNtuples)

# process.p = cms.Path(process.Stage2FileReader)
