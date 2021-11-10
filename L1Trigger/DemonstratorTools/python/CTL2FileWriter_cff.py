import FWCore.ParameterSet.Config as cms

CTL2FileWriter = cms.EDAnalyzer('CTL2FileWriter',
  puppiCands = cms.untracked.InputTag("l1ctLayer1", "Puppi"),
  inputFilename = cms.untracked.string("CTL2InputFile"),
  jets = cms.untracked.InputTag("Phase1L1TJetProducer", "UncalibratedPhase1L1TJetFromPfCandidates", ""),
  met = cms.untracked.InputTag("Phase1L1TJetProducer", "UncalibratedPhase1L1TJetFromPfCandidatesMET", ""),
  ht = cms.untracked.InputTag("Phase1L1TJetSumsProducer", "Sums", ""),
  outputFilename = cms.untracked.string("CTL2OutputFile"),
  format = cms.untracked.string("EMP")
)