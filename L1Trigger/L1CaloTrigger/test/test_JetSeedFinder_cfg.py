import FWCore.ParameterSet.Config as cms
import FWCore.Utilities.FileUtils as FileUtils # ADDED

process = cms.Process("TEST")

process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 1000

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(100) )

# fileList = FileUtils.loadListFromFile('ttbar_jets.txt')
fileList = FileUtils.loadListFromFile('ttbar_jets.txt')
readFiles = cms.untracked.vstring(*fileList)

process.source = process.source = cms.Source("PoolSource",
  fileNames = readFiles,
  # fileNames = cms.untracked.vstring(
  #  "file:debugInputs.root",
  # ),
# skipEvents = cms.untracked.uint32(3)
)

# Loads jet seeding sequence
process.load('L1Trigger.L1CaloTrigger.Phase1L1TJetSeeds_cff')

process.out = cms.OutputModule("PoolOutputModule",
  fileName = cms.untracked.string('myOutputFile_ttbar.root'),
  outputCommands = cms.untracked.vstring(
    "drop *",
    "keep *_Phase1L1TJetSeedProducer9x9trimmed_*_*",
    "keep *_l1ctLayer1_*_*"
  ),
)

process.p = cms.Path(process.Phase1L1TJetSeedsSequence )

process.e = cms.EndPath(process.out)
