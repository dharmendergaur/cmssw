import FWCore.ParameterSet.Config as cms

l1GeneratorTree = cms.EDAnalyzer(
    "L1GenTreeProducer",
    genJetToken     = cms.untracked.InputTag("ak4GenJetsNoNu"),
    jetFlavourInfosToken = cms.untracked.InputTag("slimmedGenJetsFlavourInfos"),
    genParticleToken = cms.untracked.InputTag("genParticles"),
    pileupInfoToken     = cms.untracked.InputTag("addPileupInfo"),
    genInfoToken     = cms.InputTag("generator")
)
