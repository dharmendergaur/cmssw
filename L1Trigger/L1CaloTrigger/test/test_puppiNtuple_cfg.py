import FWCore.ParameterSet.Config as cms
import FWCore.Utilities.FileUtils as FileUtils # ADDED

process = cms.Process("TEST")

process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 1000

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

sample = 'HH4B_withExtended'

fileList = FileUtils.loadListFromFile('{sample}.list'.format(sample=sample))
readFiles = cms.untracked.vstring(*fileList)

process.source = process.source = cms.Source("PoolSource",
  fileNames = readFiles,
)

process.TFileService = cms.Service(
    "TFileService",
    fileName = cms.string('puppiNtuple_{sample}.root'.format(sample=sample))
)

process.puppiNtuples = cms.EDAnalyzer('PuppiNtuples',
  inputPFCands = cms.untracked.InputTag('l1ctLayer1Extended', 'Puppi'),
)


process.p = cms.Path(process.puppiNtuples)
