import FWCore.ParameterSet.Config as cms

CTL2JetSeedsFileWriter = cms.EDAnalyzer('CTL2JetSeedWriter',
  puppiCands = cms.untracked.InputTag("l1ctLayer1", "Puppi"),
  inputFilename = cms.untracked.string("CTL2JetSeedsInputFile"),
  jetSeeds = cms.untracked.InputTag("Phase1L1TJetSeedProducer9x9trimmed", "histoJetSeeds9x9trimmed", ""),
  outputFilename = cms.untracked.string("CTL2JetSeedsOutputFile"),
  format = cms.untracked.string("EMP")
)