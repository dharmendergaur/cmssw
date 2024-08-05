import FWCore.ParameterSet.Config as cms
import FWCore.Utilities.FileUtils as FileUtils
import FWCore.ParameterSet.VarParsing as VarParsing
import math

# PART 1 : PARSE ARGUMENTS

options = VarParsing.VarParsing ('analysis')
options.register ('format',
                  'EMPv2', # default value
                  VarParsing.VarParsing.multiplicity.singleton,
                  VarParsing.VarParsing.varType.string,
                  "File format (APx, EMP or X20)")
options.parseArguments()

inputFiles = []
for filePath in options.inputFiles:
    if filePath.endswith(".root"):
        inputFiles.append(filePath)
    elif filePath.endswith("_cff.py"):
        inputFilesImport = getattr(__import__(filePath.strip(".py"),fromlist=["readFiles"]),"readFiles")
        inputFiles.extend( inputFilesImport )
    else:
        inputFiles += FileUtils.loadListFromFile(filePath)

# PART 2: SETUP MAIN CMSSW PROCESS 

process = cms.Process("GTTFileWriter")

process.load('Configuration.Geometry.GeometryExtended2026D88Reco_cff')
process.load('Configuration.Geometry.GeometryExtended2026D88_cff')
process.load('Configuration.StandardSequences.MagneticField_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:upgradePLS3', '')
process.load("FWCore.MessageLogger.MessageLogger_cfi")

process.source = cms.Source("PoolSource", fileNames = cms.untracked.vstring(inputFiles) )
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(options.maxEvents) )

process.load('SimCalorimetry.HcalTrigPrimProducers.hcaltpdigi_cff') # needed to read HCal TPs
process.load('SimCalorimetry.HGCalSimProducers.hgcalDigitizer_cfi') # needed for HGCAL_noise_fC
process.load('Configuration.StandardSequences.SimL1Emulator_cff')
process.l1t9x9HistoSeedsSCPFL1PuppiCorrectedEmulatorDC.doCorrections = True
from L1Trigger.Phase2L1ParticleFlow.l1tMHTPFProducer_cfi import l1tMHTPFProducer
process.l1t9x9HistoSeedsSCPFL1PuppiCorrectedEmulatorDCMHT = l1tMHTPFProducer.clone(jets = 'l1t9x9HistoSeedsSCPFL1PuppiCorrectedEmulatorDC')
process.extraPFStuff = cms.Task(
        process.l1tPhase2L1CaloEGammaEmulator,
        process.l1tSAMuonsGmt,
        process.l1tGTTInputProducer,
        process.l1tTrackSelectionProducer,
        process.l1tVertexFinderEmulator,
        process.L1TLayer1TaskInputsTask,
        process.L1TLayer1Task,
        process.L1TLayer2EGTask,
        process.L1TPF9x9HistoSeedJetsTaskDC,
        process.l1t9x9HistoSeedsSCPFL1PuppiCorrectedEmulatorDCMHT,
        process.L1TPF9x9HistoSeed8JetsTaskDC)

process.l1tTrackSelectionProducer.processSimulatedTracks = False

process.l1tPhase1JetSeedProducer9x9DC.phiRegions = tuple([-math.pi+i*2*math.pi/9 for i in range(10)])

process.load('L1Trigger.DemonstratorTools.CTL2FileWriter_cff')

process.CTL2FileWriter.format = cms.untracked.string(options.format)

from L1Trigger.Phase2L1ParticleFlow.l1tJetFileWriter_cfi import l1tSeededConeJetFileWriter
l1ctLayer2SCJetsProducts = cms.VPSet([cms.PSet(jets = cms.InputTag("l1t9x9HistoSeedsSCPFL1PuppiCorrectedEmulatorDC"),
                                               nJets = cms.uint32(12),
                                               mht  = cms.InputTag("l1t9x9HistoSeedsSCPFL1PuppiCorrectedEmulatorDCMHT"),
                                               nSums = cms.uint32(2)),
                                      cms.PSet(jets = cms.InputTag("l1t9x9HistoSeedsSC8PFL1PuppiCorrectedEmulatorDC"),
                                               nJets = cms.uint32(12))
                                      ])
process.l1tLayer2SeedConeJetWriter = l1tSeededConeJetFileWriter.clone(collections = l1ctLayer2SCJetsProducts)

process.MessageLogger.cerr.FwkReport.reportEvery = 1
process.Timing = cms.Service("Timing", summaryOnly = cms.untracked.bool(True))

process.p = cms.Path(process.CTL2FileWriter*process.l1tLayer2SeedConeJetWriter)
process.p.associate(process.extraPFStuff)
