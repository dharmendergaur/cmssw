import FWCore.ParameterSet.Config as cms


l1PhaseIITree = cms.EDAnalyzer("L1PhaseIITreeStep1Producer",

   egTokenBarrel = cms.InputTag("l1tEGammaClusterEmuProducer",""),  #keep as is, not produced by GCT
   tkEGTokenBarrel = cms.InputTag("l1tLayer1EG","L1TkEleEB"),
   tkEMTokenBarrel = cms.InputTag("l1tLayer1EG","L1TkEmEB"),

   egTokenHGC = cms.InputTag("l1tLayer1EG","L1EgEE"),
   tkEGTokenHGC = cms.InputTag("l1tLayer1EG","L1TkEleEE"),
   tkEMTokenHGC = cms.InputTag("l1tLayer1EG","L1TkEmEE"),

#   muonKalman = cms.InputTag("simKBmtfDigis","BMTF"), # we should remove all of these 
#   muonOverlap = cms.InputTag("simOmtfDigis","OMTF"),
#   muonEndcap = cms.InputTag("simEmtfDigis",""),
#   TkMuonToken = cms.InputTag(""),#"L1TkMuons",""),  # removing this from the run 

   #Global muons
#   muonToken = cms.untracked.InputTag("simGmtStage2Digis", ""),
#   TkGlbMuonToken = cms.InputTag(""),#L1TkGlbMuons",""), # removing this from the run 

   #GMT muons
   gmtMuonToken = cms.InputTag("l1tSAMuonsGmt", "promptSAMuons"), #we use the prompt
   gmtTkMuonToken = cms.InputTag("l1tTkMuonsGmt",""),


   scPFL1Puppi = cms.InputTag("l1tSCPFL1PuppiCorrectedEmulator", ""), #Emulator and corrected JEC; seededcone jets
   scPFL1PuppiMHT = cms.InputTag("l1tSCPFL1PuppiCorrectedEmulatorMHT", ""), #Emulator for seededcone puppiMHT

   l1pfPhase1L1TJetToken  = cms.InputTag("l1tPhase1JetCalibrator9x9" ,   "Phase1L1TJetFromPfCandidates"), #use the 9x9 case
   l1pfPhase1L1TJetMET  = cms.InputTag("l1tPhase1JetProducer9x9" ,   "UncalibratedPhase1L1TJetFromPfCandidatesMET"), #use the 9x9 case
   l1pfPhase1L1TJetSums  = cms.InputTag("l1tPhase1JetSumsProducer9x9",  "Sums"), #use the 9x9 case

   caloJetToken = cms.InputTag("l1tCaloJet","L1CaloJetCollectionBXV"),
   caloJetHTTToken = cms.InputTag("l1tCaloJetHTT","CaloJetHTT"),
   caloTauToken = cms.InputTag("l1tCaloJet","L1CaloTauCollectionBXV"),
   L1HPSPFTauToken = cms.InputTag("l1tHPSPFTauProducerPF",""),

   l1PFMet = cms.InputTag("l1tMETPFProducer",""), #emulator

   #zoPuppi = cms.InputTag(""), # does not exist anymore! 
   #l1vertextdr = cms.InputTag("VertexProducer","l1vertextdr"), #not used anymore - but kept in the loop just to be sure, not filled to ntuples
   #l1vertices = cms.InputTag("VertexProducer","l1vertices"), #not used anymore - but kept in the loop just to be sure, not filled to ntuples
   l1TkPrimaryVertex= cms.InputTag("l1tVertexFinderEmulator","l1verticesEmulation"), #we need to rename this, but these are now emulated vertices!

   L1NNTauToken = cms.InputTag("l1tNNTauProducerPuppi","L1PFTausNN"), # default collection, emulated 
   L1NNTau2vtxToken = cms.InputTag("l1tNNTauProducerPuppi2Vtx","L1PFTausNN"), # 2 vtx version 

   tkTrackerJetToken = cms.InputTag("l1tTrackJetsEmulation", "L1TrackJets"),  #these are emulated
   tkTrackerJetDisplacedToken = cms.InputTag("l1tTrackJetsExtendedEmulation", "L1TrackJetsExtended"), #emulated
	 
   tkMetToken = cms.InputTag("l1tTrackerEmuEtMiss","L1TrackerEmuEtMiss"), #emulated

   tkMhtToken = cms.InputTag("l1tTrackerEmuHTMiss","L1TrackerEmuHTMiss"), #emulated
   tkMhtDisplacedToken = cms.InputTag("l1tTrackerEmuHTMissExtended","L1TrackerEmuHTMissExtended"), #emulated

   maxL1Extra = cms.uint32(50)
)

#### Gen level tree

from L1Trigger.L1TNtuples.l1GeneratorTree_cfi  import l1GeneratorTree
genTree=l1GeneratorTree.clone()

runmenutree=cms.Path(l1PhaseIITree*genTree)




