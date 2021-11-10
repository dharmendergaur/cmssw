import FWCore.ParameterSet.Config as cms

from L1Trigger.L1CaloTrigger.Phase1L1TJetProducer_cfi import Phase1L1TJetProducer
from L1Trigger.L1CaloTrigger.Phase1L1TJetCalibrator_cfi import Phase1L1TJetCalibrator
from L1Trigger.L1CaloTrigger.Phase1L1TJetSumsProducer_cfi import Phase1L1TJetSumsProducer

Phase1L1TJetCalibratorUnsorted = Phase1L1TJetCalibrator.clone(
	  inputCollectionTag = cms.InputTag("Phase1L1TJetProducer", "UncalibratedPhase1L1TJetFromPfCandidatesUnsorted", ""),
    outputCollectionName = cms.string("Phase1L1TJetFromPfCandidates")
)

Phase1L1TJetsSequence = cms.Sequence(
  Phase1L1TJetProducer +
  Phase1L1TJetCalibrator +
  Phase1L1TJetCalibratorUnsorted +
  Phase1L1TJetSumsProducer
)
