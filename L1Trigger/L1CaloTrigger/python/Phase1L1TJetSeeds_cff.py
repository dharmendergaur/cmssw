import FWCore.ParameterSet.Config as cms

from L1Trigger.L1CaloTrigger.Phase1L1TJetSeedProducer_cfi import Phase1L1TJetSeedProducer9x9trimmed

Phase1L1TJetSeedsSequence = cms.Sequence(
  Phase1L1TJetSeedProducer9x9trimmed
)
