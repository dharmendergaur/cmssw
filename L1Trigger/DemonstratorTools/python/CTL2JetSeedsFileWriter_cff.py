import FWCore.ParameterSet.Config as cms

CTL2JetSeedsFileWriter = cms.EDAnalyzer('CTL2JetSeedWriter',
  puppiCands = cms.untracked.InputTag("l1tLayer1", "Puppi"),
  inputFilename = cms.untracked.string("CTL2JetSeedsInputFile"),
  jetSeeds = cms.untracked.InputTag("seeds9x9Trimmed", "histoJetSeeds9x9trimmed", ""),
  seededConeJets = cms.untracked.InputTag("seededcone", "", ""),
  outputFilename = cms.untracked.string("CTL2JetSeedsOutputFile"),
  format = cms.untracked.string("EMPv2")
)