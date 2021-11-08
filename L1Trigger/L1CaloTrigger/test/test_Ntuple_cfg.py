import FWCore.ParameterSet.Config as cms
from math import pi
import FWCore.Utilities.FileUtils as FileUtils
import FWCore.ParameterSet.VarParsing as VarParsing


process = cms.Process("Ntuples")

process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 100

options = VarParsing.VarParsing ('analysis')
# get and parse the command line arguments

options.register('skipEvents',
                 0,
                 VarParsing.VarParsing.multiplicity.singleton,
                 VarParsing.VarParsing.varType.int,
                 "Number of events to skip")
options.register('outFile',
                 'L1Ntuple.root',
                  VarParsing.VarParsing.multiplicity.singleton,
                 VarParsing.VarParsing.varType.string,
                 'Output file')
options.register('sample',
                 'ttbar',
                  VarParsing.VarParsing.multiplicity.singleton,
                 VarParsing.VarParsing.varType.string,
                 'Sample')
options.register('jet',
                 '9x9',
                  VarParsing.VarParsing.multiplicity.singleton,
                 VarParsing.VarParsing.varType.string,
                 'Jet')

options.parseArguments()

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(options.maxEvents)
)

fileList = FileUtils.loadListFromFile('{sample}_jets.txt'.format(sample=options.sample))
readFiles = cms.untracked.vstring(*fileList)

process.source = process.source = cms.Source("PoolSource",
  fileNames = readFiles
  # fileNames = cms.untracked.vstring( 'file:myOutputFile.root' ),
  # duplicateCheckMode = cms.untracked.string('noDuplicateCheck'),
)



process.load("L1Trigger.L1TNtuples.l1PhaseIPFJetTreeProducer_cfi")

if options.jet == '9x9':
  process.l1PhaseIPFJetTree.l1PhaseIPFJets = cms.untracked.InputTag("Phase1L1TJetCalibrator9x9", "Phase1L1TJetFromPfCandidates")
  process.l1PhaseIPFJetTree.phaseIL1PFJetSums = cms.untracked.InputTag("Phase1L1TJetSumsProducer9x9", "Sums")

if options.jet == '9x9trimmed':
  process.l1PhaseIPFJetTree.l1PhaseIPFJets = cms.untracked.InputTag("Phase1L1TJetCalibrator9x9trimmed", "Phase1L1TJetFromPfCandidates")
  process.l1PhaseIPFJetTree.phaseIL1PFJetSums = cms.untracked.InputTag("Phase1L1TJetSumsProducer9x9trimmed", "Sums")


process.TFileService = cms.Service("TFileService",
    fileName = cms.string('L1Ntuple_{sample}_{jet}.root'.format(sample = options.sample, jet = options.jet))
)

# process.p = cms.Path(process.l1PhaseIPFJetTree+process.l1PhaseIPFJetTree9x9+process.l1PhaseIPFJetTree9x9trimmed)
process.p = cms.Path(process.l1PhaseIPFJetTree)
