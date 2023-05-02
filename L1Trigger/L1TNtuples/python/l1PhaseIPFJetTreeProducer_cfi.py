import FWCore.ParameterSet.Config as cms

l1PhaseIPFJetTree = cms.EDAnalyzer("L1PhaseIPFJetTreeProducer",
   genJetToken     = cms.untracked.InputTag("ak4GenJetsNoNu"),
   jetFlavourInfosToken = cms.untracked.InputTag("slimmedGenJetsFlavourInfos"),
   genMetToken     = cms.untracked.InputTag("genMetTrue"),
   ak4L1PF = cms.untracked.InputTag("slimmedGenJetsFlavourInfos"),
   pfMetToken = cms.untracked.InputTag("l1PFMetPuppi"),
   l1PhaseIPFJets = cms.untracked.InputTag("Phase1L1TJetCalibrator", "Phase1L1TJetFromPfCandidates"),
   phaseIL1PFJetSums = cms.untracked.InputTag("Phase1L1TJetSumsProducer", "Sums"),
   phaseIL1PFMET = cms.untracked.InputTag("Phase1L1TJetProducer", "UncalibratedPhase1L1TJetFromPfCandidatesMET"),
#   ak4L1PF = cms.InputTag("L1TCorrectedPFJetProducer", "ak4PFL1PuppiCorrected"),
   cnnScore = cms.untracked.InputTag("CNNProducer", "KerasCNN"),
   maxL1Extra = cms.uint32(20)
)

runmenutree=cms.Path(l1PhaseIPFJetTree)




