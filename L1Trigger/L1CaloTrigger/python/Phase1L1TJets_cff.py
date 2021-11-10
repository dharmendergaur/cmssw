import FWCore.ParameterSet.Config as cms

from L1Trigger.L1CaloTrigger.l1tPhase1JetProducer_cfi import l1tPhase1JetProducer
from L1Trigger.L1CaloTrigger.l1tPhase1JetCalibrator_cfi import l1tPhase1JetCalibrator
from L1Trigger.L1CaloTrigger.l1tPhase1JetSumsProducer_cfi import l1tPhase1JetSumsProducer

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
