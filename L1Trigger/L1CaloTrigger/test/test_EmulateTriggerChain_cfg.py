import FWCore.ParameterSet.Config as cms
import FWCore.Utilities.FileUtils as FileUtils
from math import pi
from L1Trigger.L1CaloTrigger.Phase1L1TJetSumsProducer_cfi import Phase1L1TJetSumsProducer
from L1Trigger.L1CaloTrigger.Phase1L1TJetProducer_cfi import Phase1L1TJetProducer


process = cms.Process("TEST")

process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 1

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(1) )

# process.TFileService = cms.Service('TFileService', fileName = cms.string("CMSSWSums.root"))

fileList = FileUtils.loadListFromFile('ttbar.txt')
process.source = cms.Source("PoolSource",
  fileNames = cms.untracked.vstring(
    # "file:/hdfs/user/ec6821/L1TJets/Emulator/Processed_1117_fromScratch/runJets_TTbar_200PU_fix_with9x9trimmed/step1_Reprocess_TrackTrigger_L1_1.root"
    fileList
    ),
  # skipEvents = cms.untracked.uint32(900)
)

process.out = cms.OutputModule("PoolOutputModule",
  fileName = cms.untracked.string('ComputeJetsAndSums_TTBar_PU200_104XMTD.root'),
  outputCommands = cms.untracked.vstring(
    "drop *",
    "keep *_Phase1L1TJetProducer_*_*",
    "keep *_Phase1L1TJetSumsProducer*_*_*",
    "keep *_ak4GenJetsNoNu__*",
    "keep *_Phase1L1TJetCalibrator_*_*",
    "keep *_Phase1L1TSumsProducer_*_*",
    "keep *_Phase1PFL1TSumsProducer_*_*",
    "keep *_l1ctLayer1_Puppi_*"
  ),
)

#eta binning of the HLS FW
# caloEtaSegmentation = cms.vdouble(
#   0.0, 0.0833, 0.1666, 0.2499, 0.3332, 0.4165, 0.4998, 0.5831, 0.6664, 0.7497, 
#   0.833, 0.9163, 0.9996, 1.0829, 1.1662, 1.2495, 1.3328, 1.4161, 1.5
# )

# caloEtaSegmentation = cms.vdouble(
#   -1.5, -1.417, -1.333, -1.25, -1.167, -1.083, -1.0, -0.917, -0.833, 
#   -0.75, -0.667, -0.583, -0.5, -0.417, -0.333, -0.25, -0.167, -0.083, 
#   0.0, 0.083, 0.167, 0.25, 0.333, 0.417, 0.5, 0.583, 0.667,
#   0.75, 0.833, 0.917, 1.0, 1.083, 1.167, 1.25, 1.333, 1.417, 1.5
# )

# caloEtaSegmentation = cms.vdouble(
#   -3.0, -2.917, -2.833, -2.75, -2.667, -2.583, 
#   # -2.5, -2.417, -2.333, -2.25, -2.167, -2.083, -2.0, -1.917, -1.833, -1.75, -1.667, -1.583, 
#   -2.5, -2.417, -2.334, -2.251, -2.168, -2.085, -2.003, -1.920, -1.837, -1.754, -1.671, -1.588,
#   -1.5, -1.417, -1.333, -1.25, -1.167, -1.083, -1.0, -0.917, -0.833, -0.75, -0.667, -0.583, 
#   -0.5, -0.417, -0.333, -0.25, -0.167, -0.083, 
#   0.0, 0.083, 0.167, 0.25, 0.333, 0.417, 
#   0.5, 0.583, 0.667, 0.75, 0.833, 0.917, 
#   1.0, 1.083, 1.167, 1.25, 1.333, 1.417, 
#   1.5, 1.583, 1.667, 1.75, 1.833, 1.917, 2.0, 2.083, 2.167, 2.25, 2.333, 2.417, 
#   2.5, 2.583, 2.667, 2.75, 2.833, 2.917, 3.0)

caloEtaSegmentation = cms.vdouble(
-3 , -2.91491519955 , -2.83201206065 , -2.74910892175 , -2.66620578285 , -2.58330264395 ,
-2.5 , -2.41491519955 , -2.33201206065 , -2.24910892175 , -2.16620578285 , -2.08330264395 , -2.00039950505 , -1.91749636615 , -1.83459322725 , -1.75169008835 , -1.66878694945 , -1.58588381055 ,
-1.5 , -1.41491519955 , -1.33201206065 , -1.24910892175 , -1.16620578285 , -1.08330264395 ,
-1.0 , -0.91491519955 , -0.83201206065 , -0.74910892175 , -0.66620578285 , -0.58330264395 ,
-0.5 , -0.41491519955 , -0.33201206065 , -0.24910892175 , -0.16620578285 , -0.08330264395 ,
0 , 0.08508480045 , 0.16798793935 , 0.25089107825 , 0.33379421715 , 0.41669735605 ,
0.5 , 0.58508480045 , 0.66798793935 , 0.75089107825 , 0.83379421715 , 0.91669735605 ,
1 , 1.08508480045 , 1.16798793935 , 1.25089107825 , 1.33379421715 , 1.41669735605 ,
1.5 , 1.58508480045 , 1.66798793935 , 1.75089107825 , 1.83379421715 , 1.91669735605 , 1.99960049495 , 2.08250363385 , 2.16540677275 , 2.24830991165 , 2.33121305055 , 2.41411618945 ,
2.5 , 2.58508480045 , 2.66798793935 , 2.75089107825 , 2.83379421715 , 2.91669735605 , 3
)

# lut configuration, you can generate your own with test/generateConSinPhi.py 
# sinPhi = cms.vdouble(0.04374, 0.13087, 0.21701, 0.30149, 0.38365, 0.46289, 0.53858, 0.61015)
# cosPhi = cms.vdouble(0.99904, 0.9914, 0.97617, 0.95347, 0.92348, 0.88642, 0.84257, 0.79229)

# sinPhi = cms.vdouble(
#   -0.0353352962792, -0.122533930843, -0.208795013406, -0.293458528818, -0.375876685504, -0.455418871948, -0.531476481737, -0.603467570232, -0.670841307236, -0.733082191603, -0.789713995522, -0.840303408309, -0.884463351833, -0.921855942186, -0.952195074957, -0.975248614326, -0.990840169216, -0.998850442928, -0.999218145922, -0.99194046477, -0.977073083675, -0.954729758418, -0.925081445966, -0.888354996422, -0.844831417308, -0.794843723474, -0.738774389082, -0.677052421152, -0.610150077076, -0.538579251202, -0.462887558141, -0.383654142772, -0.301485248985, -0.217009581095, -0.130873493387, -0.0437360446299, 0.0437360446299, 0.130873493387, 0.217009581095, 0.301485248985, 0.383654142772, 0.462887558141, 0.538579251202, 0.610150077076, 0.677052421152, 0.738774389082, 0.794843723474, 0.844831417308, 0.888354996422, 0.925081445966, 0.954729758418, 0.977073083675, 0.99194046477, 0.999218145922, 0.998850442928, 0.990840169216, 0.975248614326, 0.952195074957, 0.921855942186, 0.884463351833, 0.840303408309, 0.789713995522, 0.733082191603, 0.670841307236, 0.603467570232, 0.531476481737, 0.455418871948, 0.375876685504, 0.293458528818, 0.208795013406, 0.122533930843, 0.0353352962792 )

# cosPhi = cms.vdouble( -0.999375513427, -0.992464324695, -0.977959427777, -0.955971804952, -0.926669691581, -0.890277288868, -0.847073048421, -0.797387541713, -0.741600930761, -0.680140059366, -0.613475187173, -0.542116391547, -0.466609664777, -0.387532736497, -0.305490653258, -0.2211111491, -0.135039842524, -0.0479352966351, 0.039536019772, 0.126704831606, 0.212904178348, 0.297474517214, 0.379768769555, 0.459157271892, 0.535032593708, 0.606814185113, 0.673952818851, 0.735934792636, 0.792285859677, 0.842574857312, 0.886417005995, 0.923476853383, 0.953470841004, 0.976169473869, 0.991399076421, 0.999043121392, 0.999043121392, 0.991399076421, 0.976169473869, 0.953470841004, 0.923476853383, 0.886417005995, 0.842574857312, 0.792285859677, 0.735934792636, 0.673952818851, 0.606814185113, 0.535032593708, 0.459157271892, 0.379768769555, 0.297474517214, 0.212904178348, 0.126704831606, 0.039536019772, -0.0479352966351, -0.135039842524, -0.2211111491, -0.305490653258, -0.387532736497, -0.466609664777, -0.542116391547, -0.613475187173, -0.680140059366, -0.741600930761, -0.797387541713, -0.847073048421, -0.890277288868, -0.926669691581, -0.955971804952, -0.977959427777, -0.992464324695, -0.999375513427 )

# These are the values actually used by the current HW
# Values in original LUT could not be stored precisely with current choice of ap_ufixeds
# So need to fix this

sinPhi = cms.vdouble( -0.03515625, -0.12109375, -0.20703125, -0.29296875, -0.375, -0.45703125, -0.53125, -0.6015625, -0.671875, -0.734375, -0.7890625, -0.83984375, -0.8828125, -0.921875, -0.953125, -0.9765625, -0.9921875, -1, -1, -0.9921875, -0.9765625, -0.953125, -0.92578125, -0.88671875, -0.84375, -0.79296875, -0.73828125, -0.67578125, -0.609375, -0.5390625, -0.4609375, -0.3828125, -0.30078125, -0.21875, -0.1328125, -0.04296875, 0.04296875, 0.1328125, 0.21875, 0.30078125, 0.3828125, 0.4609375, 0.5390625, 0.609375, 0.67578125, 0.73828125, 0.79296875, 0.84375, 0.88671875, 0.92578125, 0.953125, 0.9765625, 0.9921875, -1, -1, 0.9921875, 0.9765625, 0.953125, 0.921875, 0.8828125, 0.83984375, 0.7890625, 0.734375, 0.671875, 0.6015625, 0.53125, 0.45703125, 0.375, 0.29296875, 0.20703125, 0.12109375, 0.03515625
)

cosPhi = cms.vdouble( -1, -0.9921875, -0.9765625, -0.95703125, -0.92578125, -0.890625, -0.84765625, -0.796875, -0.7421875, -0.6796875, -0.61328125, -0.54296875, -0.46484375, -0.38671875, -0.3046875, -0.22265625, -0.13671875, -0.046875, 0.0390625, 0.125, 0.21484375, 0.296875, 0.37890625, 0.4609375, 0.53515625, 0.60546875, 0.67578125, 0.734375, 0.79296875, 0.84375, 0.88671875, 0.921875, 0.953125, 0.9765625, 0.9921875, -1, -1, 0.9921875, 0.9765625, 0.953125, 0.921875, 0.88671875, 0.84375, 0.79296875, 0.734375, 0.67578125, 0.60546875, 0.53515625, 0.4609375, 0.37890625, 0.296875, 0.21484375, 0.125, 0.0390625, -0.046875, -0.13671875, -0.22265625, -0.3046875, -0.38671875, -0.46484375, -0.54296875, -0.61328125, -0.6796875, -0.7421875, -0.796875, -0.84765625, -0.890625, -0.92578125, -0.95703125, -0.9765625, -0.9921875, -1 )

# sets up jet finder
process.Phase1L1TJetProducer = cms.EDProducer('Phase1L1TJetProducer',
  inputCollectionTag = cms.InputTag("l1ctLayer1", "Puppi"),
  etaBinning = caloEtaSegmentation,
  # nBinsPhi = cms.uint32(8),
  nBinsPhi = cms.uint32(72),
  phiLow = cms.double(-3.15),
  # phiUp = cms.double(0.7),
  phiUp = cms.double(3.15),
  jetIEtaSize = cms.uint32(9),
  jetIPhiSize = cms.uint32(9),
  trimmedGrid = cms.bool(True),
  seedPtThreshold = cms.double(5), # GeV
  puSubtraction = cms.bool(False),
  philsb = cms.double(0.0043633231),
  etalsb = cms.double(0.0043633231),
  outputCollectionName = cms.string("UncalibratedPhase1L1TJetFromPfCandidates"),
  vetoZeroPt = cms.bool(True),
  # pfEtaRegions = cms.vdouble( -5., -4.5, -4., -3.5, -3., -2.5, -1.5, -0.75, 0, 0.75, 1.5, 2.5, 3., 3.5, 4., 4.5, 5. ),
  etaRegions = cms.vdouble( -3, -2.5, -1.5, -1.0, -0.5, 0, 0.5, 1, 1.5, 2.5, 3 ),
  phiRegions = cms.vdouble( -3.15, -2.45, -1.75, -1.05, -0.35, 0.35, 1.05, 1.75, 2.45, 3.15 ),#, 4.2, 4.9, 5.6, 6.3 ),
  maxInputsPerRegion = cms.uint32( 18 ),
  sinPhi = sinPhi,
  cosPhi = cosPhi,
  metAbsEtaCut = cms.double(3),
  metHFAbsEtaCut = cms.double(5)
)

# sum module
process.Phase1L1TJetSumsProducer = cms.EDProducer('Phase1L1TJetSumsProducer',
  inputJetCollectionTag = cms.InputTag("Phase1L1TJetProducer", "UncalibratedPhase1L1TJetFromPfCandidatesUnsorted"),
  nBinsPhi = cms.uint32(72),
  phiLow = cms.double(-3.15),
  phiUp = cms.double(3.15),
  sinPhi = sinPhi,
  cosPhi = cosPhi,
  htPtThreshold = cms.double(30),
  htAbsEtaCut = cms.double(3),
  mhtPtThreshold = cms.double(30),
  mhtAbsEtaCut = cms.double(3),
  outputCollectionName = cms.string("Sums")
)


# process.SaveSums = cms.EDAnalyzer("SaveGenSumsAndL1Sums",
#   genMETCollectionTag = cms.InputTag("genMetTrue"), # taking pre-existing MET collection
#   l1tMETCollectionTag = cms.InputTag("Phase1L1TJetProducer", "UncalibratedPhase1L1TJetFromPfCandidatesMET"), # taking L1T MET produced by jet trigger
#   genJetCollectionTag = cms.InputTag("ak4GenJetsNoNu"), # taking pre-existing gen jet collection
#   l1tHTCollectionTag = cms.InputTag("Phase1L1TSumsProducer", "Sums"), # taking L1T HT produced by jet trigger
#   l1tMHTCollectionTag = cms.InputTag("Phase1L1TSumsProducer", "Sums"), # taking L1T MHT produced by jet trigger
# )

# process.SaveJets = cms.EDAnalyzer(
#   "StoreCandidatesToTree", 
#   candidateCollectionTag = cms.InputTag("Phase1L1TJetProducer", "UncalibratedPhase1L1TJetFromPfCandidates"),
#   treeName = cms.string("EmulatorJets"),
#   maxNumberOfCandidates = cms.uint32(100) # demonstrator can fit up to three jets
# )

# runs the jet finder and sum producer
process.p = cms.Path(
  process.Phase1L1TJetProducer + 
  process.Phase1L1TJetSumsProducer# + 
  # process.SaveSums + 
  # process.SaveJets
)

process.e = cms.EndPath(process.out)
