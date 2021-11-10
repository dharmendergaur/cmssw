import FWCore.ParameterSet.Config as cms

from L1Trigger.L1CaloTrigger.l1tPhase1JetProducer_cfi import l1tPhase1JetProducer
from L1Trigger.L1CaloTrigger.l1tPhase1JetCalibrator_9x9trimmed_cfi import l1tPhase1JetCalibrator9x9trimmed
from L1Trigger.L1CaloTrigger.l1tPhase1JetSumsProducer_cfi import l1tPhase1JetSumsProducer

l1tPhase1JetProducer9x9trimmed = l1tPhase1JetProducer.clone(
	  jetIEtaSize = 9,
	  jetIPhiSize = 9,
	  trimmedGrid = True,
	  outputCollectionName = "Uncalibratedl1tPhase1JetFromPfCandidates"
)

l1tPhase1JetCalibrator9x9trimmed.inputCollectionTag = cms.InputTag("l1tPhase1JetProducer9x9trimmed", "Uncalibratedl1tPhase1JetFromPfCandidates", "")
l1tPhase1JetCalibrator9x9trimmed.outputCollectionName = cms.string("l1tPhase1JetFromPfCandidates")

l1tPhase1JetCalibrator9x9trimmedUnsorted = l1tPhase1JetCalibrator9x9trimmed.clone(
	  inputCollectionTag = cms.InputTag("l1tPhase1JetProducer9x9trimmed", "Uncalibratedl1tPhase1JetFromPfCandidatesUnsorted", ""),
    outputCollectionName = cms.string("l1tPhase1JetFromPfCandidates")
)

l1tPhase1JetSumsProducer9x9trimmed = l1tPhase1JetSumsProducer.clone(
  inputJetCollectionTag = cms.InputTag("l1tPhase1JetCalibrator9x9trimmed", "l1tPhase1JetFromPfCandidates"),
)

l1tPhase1JetsSequence9x9trimmed = cms.Sequence(
  l1tPhase1JetProducer9x9trimmed *
  l1tPhase1JetCalibrator9x9trimmed * 
  l1tPhase1JetCalibrator9x9trimmedUnsorted *
  l1tPhase1JetSumsProducer9x9trimmed
)

L1TPFJetsPhase1Task_9x9trimmed = cms.Task(  
  l1tPhase1JetProducer9x9trimmed,
  l1tPhase1JetCalibrator9x9trimmed,
  l1tPhase1JetCalibrator9x9trimmedUnsorted,
  l1tPhase1JetSumsProducer9x9trimmed
)
