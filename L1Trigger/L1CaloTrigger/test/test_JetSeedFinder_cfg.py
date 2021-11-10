import FWCore.ParameterSet.Config as cms
import FWCore.Utilities.FileUtils as FileUtils # ADDED

process = cms.Process("TEST")

process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 1

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(100) )

fileList = FileUtils.loadListFromFile('ttbar.txt')
readFiles = cms.untracked.vstring(*fileList)

process.source = process.source = cms.Source("PoolSource",
  # fileNames = readFiles,
  fileNames = cms.untracked.vstring(
   "file:/hdfs/user/ec6821/L1TJets/LocalInputs/TTbar_200PU_Fall22_000c5e5f-78f7-44ee-95fe-7b2f2c2e2312.root",
  ),
# skipEvents = cms.untracked.uint32(3)
)

# Loads jet seeding sequence
from L1Trigger.L1CaloTrigger.Phase1L1TJetSeedProducer_cfi import l1tPhase1JetSeedProducer9x9trimmed
process.seeds9x9Trimmed = l1tPhase1JetSeedProducer9x9trimmed
# Loads seeded cone sequence
from L1Trigger.Phase2L1ParticleFlow.l1tDeregionizerProducer_cfi import l1tDeregionizerProducer
process.deregionizer = l1tDeregionizerProducer
from L1Trigger.Phase2L1ParticleFlow.l1tSeedConePFJetProducer_cfi import l1tSeedConePFJetEmulatorProducer
process.seededcone = l1tSeedConePFJetEmulatorProducer.clone(
    debug = False,
    useExternalSeeds = True,
    JetSeeds = ('seeds9x9Trimmed', 'histoJetSeeds9x9trimmed')
)
process.seededconeNoSeeds = l1tSeedConePFJetEmulatorProducer.clone(
    debug = False,
    useExternalSeeds = False,
    JetSeeds = ('seeds9x9Trimmed', 'histoJetSeeds9x9trimmed')
)

# Plain 9x9 trimmed jets
process.load("L1Trigger.L1CaloTrigger.Phase1L1TJets_9x9trimmed_cff")
process.l1tPhase1JetSumsProducer9x9trimmed.inputJetCollectionTag = ("l1tPhase1JetProducer9x9trimmed", "Uncalibratedl1tPhase1JetFromPfCandidates")
process.l1tPhase1JetSumsProducer9x9trimmed.htAbsEtaCut = 3
process.l1tPhase1JetSumsProducer9x9trimmed.mhtAbsEtaCut = 3


process.load("L1Trigger.L1TNtuples.l1PhaseIITreeStep1Producer_cfi")
process.l1PhaseIITree.l1pfPhase1L1TJetToken = ("l1tPhase1JetCalibrator9x9trimmed","Phase1L1TJetFromPfCandidates")
process.l1PhaseIITree.l1pfPhase1L1TJetSums = ("l1tPhase1JetSumsProducer9x9trimmed","Sums")
process.TFileService = cms.Service("TFileService",
    fileName = cms.string('ntuples_rerun.root')
)


process.out = cms.OutputModule("PoolOutputModule",
  fileName = cms.untracked.string('myOutputFile_ttbar.root'),
  outputCommands = cms.untracked.vstring(
    "drop *",
    "keep *_seed*_*_*",
    "keep *_ak4GenJetsNoNu_*_*",
    "keep *_Phase1L1TJetSeedProducer9x9trimmed_*_*",
    "keep *_l1tLayer1_Puppi_*",
    "keep *_l1tPhase1Jet*_*_*"
  ),
)

process.p = cms.Path( process.l1tPhase1JetsSequence9x9trimmed * process.seeds9x9Trimmed * process.deregionizer * process.seededcone * process.seededconeNoSeeds )#* process.l1PhaseIITree * process.genTree )

process.e = cms.EndPath(process.out)
