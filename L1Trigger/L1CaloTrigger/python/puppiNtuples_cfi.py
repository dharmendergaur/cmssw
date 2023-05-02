import FWCore.ParameterSet.Config as cms

puppiNtuples = cms.EDAnalyzer('PuppiNtuples',
  inputPFCands = cms.untracked.InputTag('l1ctLayer1', 'Puppi'),
)