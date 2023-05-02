import FWCore.ParameterSet.Config as cms
import FWCore.Utilities.FileUtils as FileUtils # ADDED

process = cms.Process("TEST")

process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 100

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

sample = 'HH4B'

fileList = FileUtils.loadListFromFile('{sample}.list'.format(sample=sample))
# fileList = FileUtils.loadListFromFile('SNu.list')
readFiles = cms.untracked.vstring(*fileList)

process.source = process.source = cms.Source("PoolSource",
  fileNames = readFiles,
  # fileNames = cms.untracked.vstring(
  #  "file:/hdfs/user/ec6821/L1TJets/Emulator/Processed_1230pre4_fromScratch/HH4B_200PU_1perJob/step1_Reprocess_TrackTrigger_L1_1.root",
  #  "file:/hdfs/user/ec6821/L1TJets/Emulator/Processed_1230pre4_fromScratch/HH4B_200PU_1perJob/step1_Reprocess_TrackTrigger_L1_2.root",
  # ),
  #skipEvents=cms.untracked.uint32(362)
)
# process.source.inputCommands = cms.untracked.vstring("drop *_l1pfCandidates_*_*")
# Loads 7x7 sequence
process.load('L1Trigger.L1CaloTrigger.Phase1L1TJets_cff')

# Load 9x9 sequence
process.load('L1Trigger.L1CaloTrigger.Phase1L1TJets_9x9_cff')

# Load trimmed 9x9 sequence
process.load('L1Trigger.L1CaloTrigger.Phase1L1TJets_9x9trimmed_cff')

# AK4 PF jets
process.load('L1Trigger.Phase2L1ParticleFlow.l1pfJetMet_cff')
process.l1PFJets = cms.Sequence( process.l1PFJetsTask )

process.load('L1Trigger.L1CaloTrigger.CNN_producer_cff')

process.CNNProducer.imagesFileName = "images_{sample}_CTL1Inputs.root".format(sample=sample)
process.CNNProducer.runInference = cms.bool(True)
process.CNNProducer.writeImages = cms.bool(True)
process.CNNProducer.cnnGraph = "Qcnn_graph_bigNetwork.pb"
process.CNNProducer.maxPixelValue = 64
process.CNNProducer.maxInputsPerRegion = 100
process.CNNProducer.cnnImageSizePhi = 12
process.CNNProducer.cnnImageSizeEta = 12

# process.out = cms.OutputModule("PoolOutputModule",
#   fileName = cms.untracked.string('myOutputFile_{sample}_CTL1Inputs.root'.format(sample=sample)),
#   outputCommands = cms.untracked.vstring(
#     "drop *",
#     "keep *_Phase1L1TJetProducer*_*_*",
#     "keep *_Phase1L1TJetSumsProducer*_*_*",
#     "keep *_ak4GenJetsNoNu_*_*",
#     "keep *_Phase1L1TJetCalibrator*_*_*",
#     "keep *_ak4PFL1Puppi*_*_*",
#     "keep *_l1PFMetPuppi*_*_*",
#     "keep *_genMetTrue_*_*",
#     "keep *_CNNProducer_*_*",
#     "keep *_scPFL1*_*_*",
#     "keep *_genParticles_*_*",
#     "keep *_addPileupInfo_*_*",
#     "keep *_l1ctLayer1_Puppi_*"
#   ),
# )

process.load("L1Trigger.L1TNtuples.l1PhaseIPFJetTreeProducer_cfi")
process.l1PhaseIPFJetTree.l1PhaseIPFJets = cms.untracked.InputTag("Phase1L1TJetCalibrator9x9trimmed", "Phase1L1TJetFromPfCandidates")
process.l1PhaseIPFJetTree.phaseIL1PFJetSums = cms.untracked.InputTag("Phase1L1TJetSumsProducer9x9trimmed", "Sums")

process.load("L1Trigger.L1TNtuples.L1NtupleGEN_cff")

process.puppiNtuples = cms.EDAnalyzer('PuppiNtuples',
  inputPFCands = cms.untracked.InputTag('l1ctLayer1', 'Puppi'),
)

process.TFileService = cms.Service("TFileService",
    fileName = cms.string('L1Ntuple_CNN_CTL1Inputs_{sample}.root'.format(sample=sample))
)

process.makeThings = cms.Path(process.CNNProducerSequence * process.Phase1L1TJetsSequence * process.Phase1L1TJetsSequence9x9 * process.Phase1L1TJetsSequence9x9trimmed * process.l1PFJets)
process.makeNtuples = cms.Path(process.l1GeneratorTree * process.l1PhaseIPFJetTree )# * process.puppiNtuples)
process.schedule = cms.Schedule(process.makeThings, process.makeNtuples)

