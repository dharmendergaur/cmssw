import FWCore.ParameterSet.Config as cms

CTL2FileWriter = cms.EDAnalyzer('CTL2FileWriter',
  puppiCands = cms.untracked.InputTag("l1tLayer1", "Puppi"),
  inputFilename = cms.untracked.string("CTL2InputFile"),
  jets = cms.untracked.InputTag("l1tPhase1JetProducer9x9trimmed", "Uncalibratedl1tPhase1JetFromPfCandidates"),
  met = cms.untracked.InputTag("l1tPhase1JetProducer9x9trimmed", "Uncalibratedl1tPhase1JetFromPfCandidatesMET"),
  ht = cms.untracked.InputTag("l1tPhase1JetSumsProducer9x9trimmed", "Sums"),
  outputFilename = cms.untracked.string("CTL2OutputFile"),
  format = cms.untracked.string("EMPv2")
)