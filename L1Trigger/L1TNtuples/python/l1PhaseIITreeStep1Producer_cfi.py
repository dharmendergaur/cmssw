import FWCore.ParameterSet.Config as cms

process_name = "L1"
#process_name = "HLT"

l1PhaseIITree = cms.EDAnalyzer("L1PhaseIITreeStep1Producer",

   egTokenBarrel = cms.InputTag("l1tEGammaClusterEmuProducer","",process_name),  #keep as is, not produced by GCT
   tkEGTokenBarrel = cms.InputTag("l1tLayer1EG","L1TkEleEB",process_name),
   tkEMTokenBarrel = cms.InputTag("l1tLayer1EG","L1TkEmEB",process_name),

   egTokenHGC = cms.InputTag("l1tLayer1EG","L1EgEE",process_name),
   tkEGTokenHGC = cms.InputTag("l1tLayer1EG","L1TkEleEE",process_name),
   tkEMTokenHGC = cms.InputTag("l1tLayer1EG","L1TkEmEE",process_name),

                               #   muonKalman = cms.InputTag("simKBmtfDigis","BMTF",process_name), # we should remove all of these
                               #   muonOverlap = cms.InputTag("simOmtfDigis","OMTF",process_name),
                               #   muonEndcap = cms.InputTag("simEmtfDigis","",process_name),
                               #   TkMuonToken = cms.InputTag(""),#"L1TkMuons",process_name,""),  # removing this from the run

   #Global muons
                               #   muonToken = cms.untracked.InputTag("simGmtStage2Digis", "",process_name),
                               #   TkGlbMuonToken = cms.InputTag(""),#L1TkGlbMuons",",process_name"), # removing this from the run

   #GMT muons
   gmtMuonToken = cms.InputTag("l1tSAMuonsGmt", "promptSAMuons",process_name), #we use the prompt
   gmtTkMuonToken = cms.InputTag("l1tTkMuonsGmt","",process_name),


   scPFL1Puppi = cms.InputTag("l1tSCPFL1PuppiCorrectedEmulator", ""), #Emulator and corrected JEC; seededcone jets
   scPFL1PuppiMHT = cms.InputTag("l1tSCPFL1PuppiCorrectedEmulatorMHT", ""), #Emulator for seededcone puppiMHT
   scPFL1PuppiExtended = cms.InputTag("l1tSCPFL1PuppiExtendedCorrectedEmulator", ""), #Emulator and corrected JEC; seededcone jets from extended tracks/puppi
   scBjetNN = cms.InputTag("l1tBJetProducerPuppiCorrectedEmulator", "L1PFBJets"), #B Jet NN ID for seededcone jets from extended tracks/puppi

   l1pfPhase1L1TJetToken  = cms.InputTag("l1tPhase1JetCalibrator9x9" ,   "Phase1L1TJetFromPfCandidates",process_name), #use the 9x9 case
   l1pfPhase1L1TJetMET  = cms.InputTag("l1tPhase1JetProducer9x9" ,   "UncalibratedPhase1L1TJetFromPfCandidatesMET",process_name), #use the 9x9 case
   l1pfPhase1L1TJetSums  = cms.InputTag("l1tPhase1JetSumsProducer9x9",  "Sums",process_name), #use the 9x9 case

   caloJetToken = cms.InputTag("l1tCaloJet","L1CaloJetCollectionBXV",process_name),
   caloJetHTTToken = cms.InputTag("l1tCaloJetHTT","CaloJetHTT",process_name),
   caloTauToken = cms.InputTag("l1tCaloJet","L1CaloTauCollectionBXV",process_name),
   L1HPSPFTauToken = cms.InputTag("l1tHPSPFTauProducerPF","",process_name),

   l1PFMet = cms.InputTag("l1tMETPFProducer","",process_name), #emulator

   #zoPuppi = cms.InputTag(""), # does not exist anymore!
                               #l1vertextdr = cms.InputTag("VertexProducer",process_name,"l1vertextdr"), #not used anymore - but kept in the loop just to be sure, not filled to ntuples
                               #l1vertices = cms.InputTag("VertexProducer","l1vertices",process_name), #not used anymore - but kept in the loop just to be sure, not filled to ntuples
   l1TkPrimaryVertex= cms.InputTag("l1tVertexFinderEmulator","l1verticesEmulation",process_name), #we need to rename this, but these are now emulated vertices!

   L1NNTauToken = cms.InputTag("l1tNNTauProducerPuppi","L1PFTausNN",process_name), # default collection, emulated
   L1NNTau2vtxToken = cms.InputTag("l1tNNTauProducerPuppi2Vtx","L1PFTausNN",process_name), # 2 vtx version

   tkTrackerJetToken = cms.InputTag("l1tTrackJetsEmulation", "L1TrackJets",process_name),  #these are emulated
   tkTrackerJetDisplacedToken = cms.InputTag("l1tTrackJetsExtendedEmulation", "L1TrackJetsExtended",process_name), #emulated

   tkMetToken = cms.InputTag("l1tTrackerEmuEtMiss","L1TrackerEmuEtMiss",process_name), #emulated

   tkMhtToken = cms.InputTag("l1tTrackerEmuHTMiss","L1TrackerEmuHTMiss",process_name), #emulated
   tkMhtDisplacedToken = cms.InputTag("l1tTrackerEmuHTMissExtended","L1TrackerEmuHTMissExtended",process_name), #emulated

   maxL1Extra = cms.uint32(50)
)

#### Gen level tree

from L1Trigger.L1TNtuples.l1GeneratorTree_cfi  import l1GeneratorTree
genTree=l1GeneratorTree.clone()

runmenutree=cms.Path(l1PhaseIITree*genTree)
