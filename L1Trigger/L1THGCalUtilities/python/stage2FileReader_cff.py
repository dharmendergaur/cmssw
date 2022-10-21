import FWCore.ParameterSet.Config as cms

Stage2FileReader = cms.EDProducer('Stage2FileReader',
  files = cms.vstring("HGCS2OutputToL1TFile_Sector0_0.txt"),
  format = cms.untracked.string("EMPv2")
)
