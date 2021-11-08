import FWCore.ParameterSet.Config as cms
import FWCore.Utilities.FileUtils as FileUtils # ADDED

process = cms.Process("TEST")

process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 1

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(10) )

fileList = FileUtils.loadListFromFile('ttbar.list')
readFiles = cms.untracked.vstring(*fileList)

process.source = process.source = cms.Source("PoolSource",
  # fileNames = readFiles,
  fileNames = cms.untracked.vstring(
   "file:debugInputs.root",
  ),
# skipEvents = cms.untracked.uint32(1)
)

# Loads 7x7 sequence
process.load('L1Trigger.L1CaloTrigger.Phase1L1TJets_cff')

# Load 9x9 sequence
process.load('L1Trigger.L1CaloTrigger.Phase1L1TJets_9x9_cff')

# Load trimmed 9x9 sequence
process.load('L1Trigger.L1CaloTrigger.Phase1L1TJets_9x9trimmed_cff')

# AK4 PF jets
process.load('L1Trigger.Phase2L1ParticleFlow.l1pfJetMet_cff')
process.ak4PFL1Puppi.src = cms.InputTag("l1ctLayer1","Puppi")
process.ak4PFL1PuppiCorrected.correctorFile = cms.string('L1Trigger/Phase2L1ParticleFlow/data/jecs/jecs.PU200_110X.root')
process.l1PFJets = cms.Sequence( process.ak4PFL1Puppi + process.ak4PFL1PuppiCorrected )


process.out = cms.OutputModule("PoolOutputModule",
  fileName = cms.untracked.string('myOutputFile.root'),
  outputCommands = cms.untracked.vstring(
    "drop *",
    "keep *_Phase1L1TJetProducer*_*_*",
    "keep *_Phase1L1TJetSumsProducer*_*_*",
    "keep *_ak4GenJetsNoNu_*_*",
    "keep *_Phase1L1TJetCalibrator*_*_*",
    "keep *_ak4PFL1Puppi*_*_*",
    "keep *_l1PFMetPuppi*_*_*",
    "keep *_genMetTrue_*_*"
  ),
)

process.p = cms.Path(process.Phase1L1TJetsSequence * process.Phase1L1TJetsSequence9x9 * process.Phase1L1TJetsSequence9x9trimmed * process.l1PFJets * process.l1PFMetPuppi )

process.e = cms.EndPath(process.out)
