import FWCore.ParameterSet.Config as cms

from L1Trigger.L1CaloTrigger.l1tPhase1JetProducer_cfi import l1tPhase1JetProducer
from L1Trigger.L1CaloTrigger.l1tPhase1JetCalibrator9_cfi import l1tPhase1JetCalibrator9x9
from L1Trigger.L1CaloTrigger.l1tPhase1JetSumsProducer_cfi import l1tPhase1JetSumsProducer

l1tPhase1JetProducer9x9 = l1tPhase1JetProducer.clone(
	  jetIEtaSize = 9,
	  jetIPhiSize = 9,
	  outputCollectionName = "Uncalibratedl1tPhase1JetFromPfCandidates"
)

l1tPhase1JetCalibrator9x9.inputCollectionTag = cms.InputTag("l1tPhase1JetProducer9x9", "Uncalibratedl1tPhase1JetFromPfCandidates", "")
l1tPhase1JetCalibrator9x9.outputCollectionName = cms.string("l1tPhase1JetFromPfCandidates")

l1tPhase1JetSumsProducer9x9 = l1tPhase1JetSumsProducer.clone(
  inputJetCollectionTag = cms.InputTag("l1tPhase1JetCalibrator9x9", "l1tPhase1JetFromPfCandidates"),
)

l1tPhase1JetCalibrator9x9Unsorted = l1tPhase1JetCalibrator9x9.clone(
	  inputCollectionTag = cms.InputTag("l1tPhase1JetProducer9x9", "Uncalibratedl1tPhase1JetFromPfCandidatesUnsorted", ""),
    outputCollectionName = cms.string("l1tPhase1JetFromPfCandidates")
)


l1tPhase1JetsSequence9x9 = cms.Sequence(
  l1tPhase1JetProducer9x9 *
  l1tPhase1JetCalibrator9x9 * 
  l1tPhase1JetCalibrator9x9Unsorted *
  l1tPhase1JetSumsProducer9x9
)

L1TPFJetsPhase1Task_9x9 = cms.Task( 
    l1tPhase1JetProducer9x9,
    l1tPhase1JetCalibrator9x9,
    l1tPhase1JetCalibrator9x9Unsorted,
    l1tPhase1JetSumsProducer9x9
  )
